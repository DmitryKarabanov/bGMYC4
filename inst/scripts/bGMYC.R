#═══════════════════════════════════════════════════════════════
#bGMYC4 v4.1.0: Полный интерактивный пайплайн
#(параметры -> тест -> мульти/сингл -> визуализация с posterior-colored деревом)
#═══════════════════════════════════════════════════════════════

library(ape)
library(future)
library(future.apply)
library(mcmcse)
library(plotly)
library(base64enc)
library(htmlwidgets)
library(treeio)
library(dplyr)
library(ggtree)
library(bGMYC4)

cat("✅ Все критичные зависимости загружены\n")
cat("=== 🌿 bGMYC4 Interactive Analysis ===\n\n")

# 1. ЗАГРУЗКА ДЕРЕВЬЕВ  
consensus_path <- readline(prompt = "📥 Путь к консенсусному дереву (tree): ")
posterior_path <- readline(prompt = "📥 Путь к набору деревьев (trees): ")

if (!file.exists(consensus_path)) stop("❌ Консенсусное дерево не найдено!")
if (!file.exists(posterior_path)) stop("❌ Файл с деревьями не найден!")

DELIM_DIR <- dirname(consensus_path)
html_path <- file.path(DELIM_DIR, "bGMYC_interactive_heatmap.html")
cat(sprintf("💾 Результаты будут сохранены в: %s\n", DELIM_DIR))

# Используем treeio::read.beast() для извлечения posterior из Nexus-аннотаций
tree_beast <- read.beast(consensus_path)
tree <- tree_beast@phylo  # Извлекаем phylo-объект для bGMYC

# Для posterior_trees используем обычный read.nexus
all_trees <- read.nexus(posterior_path)
class(all_trees) <- "multiPhylo"

cat(sprintf("✅ Загружено: 1 консенсус + %d деревьев в апостериоре\n", length(all_trees)))

#  2. ВЫБОР РЕЖИМА АНАЛИЗА 
cat("\n📊 Выберите режим анализа:\n")
cat("   1. Только консенсусное дерево (быстрый анализ, без учёта топологической неопределённости)\n")
cat("   2. Консенсус + множество деревьев (учёт филогенетической неопределённости)\n")

mode_str <- readline("   Введите 1 или 2 (по умолчанию 2): ")
analysis_mode <- ifelse(nchar(trimws(mode_str)) == 0 || is.na(as.integer(mode_str)), 2, as.integer(mode_str))

if (!analysis_mode %in% c(1, 2)) {
  cat("⚠️ Неверный выбор. Установлен режим 2.\n")
  analysis_mode <- 2
}

# 3. ВЫБОРКА ДЕРЕВЬЕВ (ТОЛЬКО ДЛЯ РЕЖИМА 2)
if (analysis_mode == 2) {
  n_total <- length(all_trees)
  burnin_idx <- floor(n_total * 0.10)
  
  n_str <- readline(prompt = sprintf("\n🎲 Сколько деревьев выбрать? (по умолчанию 10, доступно: %d): ", n_total - burnin_idx))
  n_sample <- ifelse(nchar(n_str) == 0 || is.na(as.numeric(n_str)), 10, as.integer(n_str))
  
  if (n_sample > (n_total - burnin_idx)) n_sample <- n_total - burnin_idx
  
  set.seed(42)
  trees_sample <- all_trees[sample((burnin_idx + 1):n_total, n_sample)]
  class(trees_sample) <- "multiPhylo"
  
  cat(sprintf("✅ Выбрано %d случайных деревьев (индексы %d–%d)\n", n_sample, burnin_idx + 1, n_total))
} else {
  cat("\n📦 Режим 1 выбран: анализ только на консенсусном дереве.\n")
  trees_sample <- NULL
}

# 4. ПРОВЕРКА УЛЬТРАМЕТРИЧНОСТИ
fix_ultrametric <- function(tr) {
  if (!is.ultrametric(tr)) {
    tr$edge.length <- round(tr$edge.length, 8)
    if (!is.ultrametric(tr)) stop("❌ Дерево остаётся неультраметричным после округления ветвей.\n")
  }
  tr
}

cat("📐 Проверка ультраметричности...\n")
tree <- fix_ultrametric(tree)

if (analysis_mode == 2) {
  trees_sample <- lapply(trees_sample, fix_ultrametric)
  class(trees_sample) <- "multiPhylo"
}

cat("✅ Деревья ультраметричны.\n")

# 5. АУТГРУППА (С ПАТТЕРНАМИ)
og_prompt <- readline(prompt = "\n🌳 Удалять аутгруппы? (y/n, по умолчанию n): ")

if (tolower(og_prompt) == "y") {
  og_str <- readline(prompt = "   Введите имена аутгрупп (через запятую): ")
  og_patterns <- trimws(unlist(strsplit(og_str, "[,;]+")))
  og_patterns <- og_patterns[nchar(og_patterns) > 0]
  
  if (length(og_patterns) > 0) {
    matching_tips <- c()
    unmatched_patterns <- c()
    
    for(p in og_patterns) {
      pattern <- paste0("(^|[^A-Za-z0-9])", p, "($|[^A-Za-z0-9])")
      found <- tree$tip.label[grepl(pattern, tree$tip.label, ignore.case = TRUE, perl = TRUE)]
      
      if (length(found) > 0) {
        matching_tips <- c(matching_tips, found)
      } else {
        unmatched_patterns <- c(unmatched_patterns, p)
      }
    }
    
    matching_tips <- unique(matching_tips)
    
    if (length(unmatched_patterns) > 0) {
      cat(sprintf("⚠️  Не найдены в дереве: %s\n", paste(unmatched_patterns, collapse = ", ")))
    }
    
    if (length(matching_tips) > 0) {
      tree <- drop.tip(tree, matching_tips)
      
      if (exists("analysis_mode") && exists("trees_sample")) {
        if (analysis_mode == 2) {
          trees_sample <- lapply(trees_sample, drop.tip, tip = matching_tips)
          class(trees_sample) <- "multiPhylo"
        }
      }
      
      cat(sprintf("✅ Удалено %d таксонов:\n", length(matching_tips)))
      if (length(matching_tips) > 10) {
        cat(sprintf("   %s ... и еще %d\n", 
                    paste(head(matching_tips, 10), collapse = ", "), 
                    length(matching_tips) - 10))
      } else {
        cat(sprintf("   %s\n", paste(matching_tips, collapse = ", ")))
      }
    } else {
      cat("⚠️  Ни один таксон не подходит под указанные имена. Пропускаю.\n")
    }
  }
}

ntips <- length(tree$tip.label)
cat(sprintf("📊 Итоговый размер: %d таксонов\n\n", ntips))

# 6. БЛОК ЗАДАНИЯ ПАРАМЕТРОВ 
input_scalar <- function(label, default_val) {
  val <- readline(sprintf("  %s (по умолчанию: %s): ", label, default_val))
  if (nchar(trimws(val)) == 0) return(default_val)
  num <- suppressWarnings(as.numeric(val))
  if (is.na(num)) {
    cat("⚠️  Введено не число. Используется значение по умолчанию.\n")
    return(default_val)
  }
  return(num)
}

input_vector <- function(label, default_vec) {
  val <- readline(sprintf("  %s (через запятую, по умолчанию: %s): ", label, paste(default_vec, collapse = ", ")))
  if (nchar(trimws(val)) == 0) return(default_vec)
  nums <- suppressWarnings(as.numeric(unlist(strsplit(val, "[,;\\s]+"))))
  nums <- nums[!is.na(nums)]
  if (length(nums) != length(default_vec)) {
    cat(sprintf("⚠️  Ожидается %d числа. Используется значение по умолчанию.\n", length(default_vec)))
    return(default_vec)
  }
  return(nums)
}

cat("📝 ЗАДАНИЕ ПАРАМЕТРОВ МОДЕЛИ:\n")

params <- list()
params$mcmc      <- input_scalar("mcmc", 10000)
params$burnin    <- input_scalar("burnin", 1000)
params$thinning  <- input_scalar("thinning", 10)
params$py1       <- input_scalar("py1", 0)
params$py2       <- input_scalar("py2", 1.5)
params$pc1       <- input_scalar("pc1", 0)
params$pc2       <- input_scalar("pc2", 2)
params$t1        <- input_scalar("t1", 2)
params$t2        <- input_scalar("t2", min(ntips - 1, 100))
params$scale     <- input_vector("scale", c(20, 10, 5))
params$start     <- input_vector("start", c(1, 1, floor((params$t1 + min(ntips-1, 50))/2)))

params$t2 <- max(params$t1 + 2, min(params$t2, ntips - 1))
params$start[3] <- max(params$t1 + 1, min(params$start[3], params$t2 - 1))

cat(sprintf("🔒 Параметры зафиксированы: t ∈ [%d, %d] | start[3] = %d\n\n", params$t1, params$t2, params$start[3]))

# 7. ИНТЕРАКТИВНЫЙ ЦИКЛ ТЕСТА SINGLEPHY
repeat {
  cat("🔄 Запуск bgmyc.singlephy на консенсусном дереве...\n")
  
  res_single <- bgmyc.singlephy(
    phylo = tree,
    mcmc = params$mcmc, burnin = params$burnin, thinning = params$thinning,
    py1 = params$py1, py2 = params$py2, pc1 = params$pc1, pc2 = params$pc2,
    t1 = params$t1, t2 = params$t2, scale = params$scale, start = params$start
  )
  
  plot(res_single)
  cat("👀 Закройте окно графиков для анализа сходимости...\n")
  Sys.sleep(1); flush.console()
  
# 🟢 ДИНАМИЧЕСКАЯ АНАЛИТИКА СХОДИМОСТИ 
  ar <- res_single$accept
  cat(sprintf("\n📊 Acceptance rates: py=%.3f | pc=%.3f | th=%.3f\n", ar[1], ar[2], ar[3]))
  
  if (requireNamespace("mcmcse", quietly = TRUE)) {
    ess_vals <- sapply(1:4, function(col) round(mcmcse::ess(res_single$par[, col])))
    cat(sprintf("🔢 ESS: py=%d | pc=%d | th=%d | logL=%d [Оптимум: >200]\n",
                ess_vals[1], ess_vals[2], ess_vals[3], ess_vals[4]))
  } else {
    cat("⚠️  ESS: пакет 'mcmcse' не установлен. Выполните install.packages('mcmcse')\n")
    ess_vals <- NULL
  }
  
  recs <- character(0)
  param_names <- c("py", "pc", "th")
  
  for (i in 1:3) {
    if (ar[i] > 0.55) {
      new_val <- round(params$scale[i] * 1.5)
      if (new_val == params$scale[i]) new_val <- params$scale[i] + (if(i < 3) 5 else 2)
      recs <- c(recs, sprintf("• %s: слишком высокая (%.2f) → увеличьте scale[%d] с %g до ~%d",
                              param_names[i], ar[i], i, params$scale[i], new_val))
    } else if (ar[i] < 0.15) {
      new_val <- max(if(i < 3) 2 else 1, round(params$scale[i] * 0.5))
      recs <- c(recs, sprintf("• %s: слишком низкая (%.2f) → уменьшите scale[%d] с %g до ~%d",
                              param_names[i], ar[i], i, params$scale[i], new_val))
    }
  }
  
  if (!is.null(ess_vals) && any(ess_vals < 200)) {
    recs <- c(recs, "• ↑ mcmc или ↓ thinning (ESS < 200: цепь требует больше независимых шагов)")
  }
  
  if (length(recs) > 0) {
    cat("\n💡 РЕКОМЕНДАЦИИ (целевой диапазон 0.20–0.40):\n")
    for (r in recs) cat(sprintf("   %s\n", r))
    cat("   💡 При повторном запуске теста подставьте предложенные значения в поле scale.\n")
  } else {
    cat("✅ Сходимость цепи стабильна. Параметры оптимальны.\n")
  }
  

  next_prompt <- if (analysis_mode == 1) {
    "⏭️  Перейти к финальному выводу (y) или изменить параметры (n)? [y/n]: "
  } else {
    "⏭️  Перейти к multiphylo (y) или изменить параметры (n)? [y/n]: "
  }
  
  choice <- readline(prompt = sprintf("\n%s", next_prompt))
  if (tolower(choice) != "n") break
  
  cat("\n📝 Обновите параметры (Enter = оставить текущее):\n")
  
  params$mcmc      <- input_scalar("mcmc", params$mcmc)
  params$burnin    <- input_scalar("burnin", params$burnin)
  params$thinning  <- input_scalar("thinning", params$thinning)
  params$py1       <- input_scalar("py1", params$py1)
  params$py2       <- input_scalar("py2", params$py2)
  params$pc1       <- input_scalar("pc1", params$pc1)
  params$pc2       <- input_scalar("pc2", params$pc2)
  params$t1        <- input_scalar("t1", params$t1)
  params$t2        <- input_scalar("t2", params$t2)
  params$scale     <- input_vector("scale", params$scale)
  params$start     <- input_vector("start", params$start)
  
  if (length(params$scale) != 3) params$scale <- c(20, 10, 5)
  if (length(params$start) != 3) params$start <- c(1, 0.5, floor((params$t1 + params$t2)/2))
  
  params$t2 <- max(params$t1 + 2, min(params$t2, ntips - 1))
  params$start[3] <- max(params$t1 + 1, min(params$start[3], params$t2 - 1))
  
  cat(sprintf("🔒 Параметры обновлены: t ∈ [%d, %d] | start[3] = %d\n\n", params$t1, params$t2, params$start[3]))
}

# 8. ПОДГОТОВКА ФИНАЛЬНЫХ РЕЗУЛЬТАТОВ
if (analysis_mode == 1) {
  cat("\n🌲 Анализ на одном дереве завершён. Формирую выводы...\n")
  final_res <- list(res_single)
  class(final_res) <- "multibgmyc"
} else {
  cat("\n⚡ Запуск bgmyc.multiphylo на выбранных деревьях...\n")
  
  n_physical <- parallel::detectCores(logical = FALSE)
  n_workers <- min(n_physical - 1, n_sample)
  if (n_workers < 1) n_workers <- 1
  
  cat(sprintf("   Воркеры: %d (физических ядер: %d)\n", n_workers, n_physical))
  
  plan(multisession, workers = n_workers)
  
  final_res <- future_lapply(seq_along(trees_sample), function(i) {
    bgmyc.singlephy(
      phylo = trees_sample[[i]],
      mcmc = params$mcmc, burnin = params$burnin, thinning = params$thinning,
      py1 = params$py1, py2 = params$py2, pc1 = params$pc1, pc2 = params$pc2,
      t1 = params$t1, t2 = params$t2, scale = params$scale, start = params$start
    )
  }, future.seed = TRUE)
  
  class(final_res) <- "multibgmyc"
  cat("✅ Все деревья успешно обработаны.\n")
  
  if (requireNamespace("mcmcse", quietly = TRUE) && length(final_res) > 1) {
    cat("\n📐 Расчёт Gelman-Rubin (R̂) across tree chains...\n")
    
    chains_list <- lapply(final_res, function(res) res$par[, 3])
    names(chains_list) <- paste0("Tree_", seq_along(final_res))
    
    gr_result <- tryCatch(mcmcse::gelman(chains_list), error = function(e) NULL)
    
    if (!is.null(gr_result)) {
      rhat <- round(gr_result$Rhat, 3)
      cat(sprintf("🔍 Gelman-Rubin R̂ (threshold): %.3f [Оптимум: < 1.05]\n", rhat))
      
      if (rhat > 1.05) {
        cat("⚠️  R̂ > 1.05: цепи показывают расхождение. Увеличьте mcmc/burnin или проверьте топологии деревьев.\n")
      } else {
        cat("✅ Цепи сошлись стабильно across posterior trees.\n")
      }
    }
  }
}

# 9. ВИЗУАЛИЗАЦИЯ И АВТОМАТИЧЕСКИЙ ЭКСПОРТ CSV 
cat("\n🌡️ Построение интерактивной карты вероятностей (CoMa-style)...\n")
probmat <- spec.probmat(final_res)

# Получаем порядок кончиков из дерева
p_tree <- suppressWarnings(ggtree(tree, layout = "rectangular"))
tips_data <- p_tree$data %>% dplyr::filter(isTip) %>% dplyr::arrange(y)
tip_order <- tips_data$label
n_tips <- length(tip_order)

# Перестраиваем матрицу под порядок дерева
if (!all(rownames(probmat) == tip_order) || !all(colnames(probmat) == tip_order)) {
  probmat <- probmat[tip_order, tip_order]
}

#  КАСТОМНАЯ ВИЗУАЛИЗАЦИЯ (дерево + матрица) 

cat("\n🎨 Создание кастомной визуализации с деревом (1:1)...\n")

# Если tree_beast не был загружен ранее, читаем его через treeio
if (!exists("tree_beast") || !inherits(tree_beast, "treedata")) {
  tree_beast <- treeio::read.beast(consensus_path)
}

# Получаем порядок кончиков из дерева
if(!is.null(tree$edge.length)) tree$edge.length[tree$edge.length < 0] <- 0
options(ignore.negative.edge = TRUE)
p_tree_raw <- suppressWarnings(ggtree(tree, layout = "rectangular"))
tree_data_raw <- p_tree_raw$data
tips_data_raw <- tree_data_raw %>% dplyr::filter(isTip) %>% dplyr::arrange(y)
tip_order <- tips_data_raw$label
n_tips_tree <- length(tip_order)

cat(sprintf("   Таксонов в дереве: %d\n", n_tips_tree))
cat(sprintf("   Таксонов в probmat: %d\n", nrow(probmat)))

# Находим пересечение таксонов между деревом и матрицей
common_taxa <- intersect(tip_order, rownames(probmat))
cat(sprintf("   Общие таксоны: %d\n", length(common_taxa)))

if (length(common_taxa) == 0) {
  stop("❌ Нет общих таксонов между деревом и матрицей вероятностей!")
}

# Фильтруем дерево только по общим таксонам, сохраняя аннотации через treeio
tips_to_drop <- setdiff(tree_beast@phylo$tip.label, common_taxa)
if(length(tips_to_drop) > 0) {
  tree_beast_filtered <- treeio::drop.tip(tree_beast, tips_to_drop)
} else {
  tree_beast_filtered <- tree_beast
}
tree_filtered <- tree_beast_filtered@phylo

p_tree <- suppressWarnings(ggtree(tree_beast_filtered, layout = "rectangular"))
tree_data <- p_tree$data
tips_data <- tree_data %>% dplyr::filter(isTip) %>% dplyr::arrange(y)
tip_order <- tips_data$label
n_tips <- length(tip_order)

# Фильтруем и переупорядочиваем матрицу
probmat <- probmat[tip_order, tip_order]

cat(sprintf("   ✅ Матрица %d×%d синхронизирована с деревом\n", n_tips, n_tips))

# Извлекаем posterior из аннотаций tree_data
post_col <- intersect(c("posterior", "prob", "Posterior", "PROB"), colnames(tree_data))
if(length(post_col) > 0) {
  post_col <- post_col[1]
  node_posterior <- tree_data[[post_col]]
  names(node_posterior) <- tree_data$node
} else {
  # Фоллбэк: парсим из label (если вдруг read.nexus оставил их там)
  node_posterior <- sapply(tree_data$label, function(lbl) {
    if(is.na(lbl) || lbl == "") return(NA_real_)
    m <- regmatches(lbl, regexpr("(posterior|prob)\\s*=\\s*([0-9\\.eE\\-]+)", lbl, ignore.case = TRUE, perl = TRUE))
    if(length(m) == 0 || m == "") return(NA_real_)
    num_str <- sub("^(posterior|prob)\\s*=\\s*", "", m, ignore.case = TRUE, perl = TRUE)
    val <- as.numeric(num_str)
    if(!is.na(val) && val > 1) val <- val / 100
    return(val)
  })
  names(node_posterior) <- tree_data$node
}


# ДЕРЕВО (левая часть)

edges <- tree_data %>% dplyr::filter(!is.na(parent))

fig_tree <- plotly::plot_ly()

# Функция цвета: плавный градиент через несколько контрольных точек
# Красный (0) -> Оранжевый (0.25) -> Жёлтый (0.5) -> Светло-зелёный (0.75) -> Зелёный (1.0)
get_pp_color <- function(pp) {
  if(is.na(pp)) return("#CCCCCC") # Серый для NA
  pp <- max(0, min(1, pp))
  
  # Определяем контрольные точки градиента
  colors <- list(
    c(0.0, 1.0, 0.0, 0.0),    # Красный
    c(0.25, 1.0, 0.5, 0.0),   # Оранжевый
    c(0.5, 1.0, 1.0, 0.0),    # Жёлтый
    c(0.75, 0.5, 1.0, 0.0),   # Светло-зелёный
    c(1.0, 0.0, 0.7, 0.0)     # Зелёный
  )
  
  # Находим два ближайших контрольных точки
  if(pp <= 0) return(sprintf("#%02X%02X%02X", 255, 0, 0))
  if(pp >= 1) return(sprintf("#%02X%02X%02X", 0, 179, 0))
  
  # Интерполяция между контрольными точками
  for(i in 1:(length(colors)-1)) {
    if(pp >= colors[[i]][1] && pp <= colors[[i+1]][1]) {
      # Линейная интерполяция между двумя точками
      t <- (pp - colors[[i]][1]) / (colors[[i+1]][1] - colors[[i]][1])
      r <- colors[[i]][2] + t * (colors[[i+1]][2] - colors[[i]][2])
      g <- colors[[i]][3] + t * (colors[[i+1]][3] - colors[[i]][3])
      b <- colors[[i]][4] + t * (colors[[i+1]][4] - colors[[i]][4])
      return(sprintf("#%02X%02X%02X", round(r*255), round(g*255), round(b*255)))
    }
  }
  
  return("#CCCCCC")
}

for(i in 1:nrow(edges)) {
  child <- edges[i, ]
  parent <- tree_data %>% dplyr::filter(node == child$parent)
  if(nrow(parent) == 0) next
  
  # Определяем цвет ветви по posterior дочернего узла
  child_node <- as.character(child$node)
  pp <- node_posterior[child_node]
  branch_color <- get_pp_color(pp)
  
  # Hover текст
  pp_text <- ifelse(is.na(pp), "N/A", sprintf("%.2f", pp))
  hover_txt <- paste0("<b>Node:</b> ", child$node, "<br>",
                      "<b>Posterior:</b> ", pp_text)
  
  # Ветви (ширина 5)
  fig_tree <- fig_tree %>% plotly::add_segments(
    x = parent$x, xend = child$x, y = child$y, yend = child$y,
    line = list(color = branch_color, width = 5), hovertext = hover_txt, hoverinfo = "text", showlegend = FALSE)
  fig_tree <- fig_tree %>% plotly::add_segments(
    x = parent$x, xend = parent$x, y = parent$y, yend = child$y,
    line = list(color = branch_color, width = 5), hovertext = hover_txt, hoverinfo = "text", showlegend = FALSE)
}

# Подписи таксонов
max_x <- max(tree_data$x, na.rm = TRUE)
label_offset <- max_x * 0.05 

# Разреженные подписи 
label_step <- max(1, ceiling(n_tips / 40))
tips_labeled <- tips_data %>%
  dplyr::mutate(row_num = dplyr::row_number()) %>%
  dplyr::filter(row_num %% label_step == 1)

fig_tree <- fig_tree %>% plotly::add_text(
  data = tips_labeled,
  x = ~x + label_offset, y = ~y, text = ~label,
  textposition = "middle left",
  textfont = list(size = 16, family = "monospace", color = "#333333"),
  showlegend = FALSE, hoverinfo = "skip")

# Настройка осей дерева (autorange = "reversed" чтобы y=1 был сверху, как в ggtree)
fig_tree <- fig_tree %>% plotly::layout(
  xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "",
               range = c(-0.02, max_x + label_offset + max_x*0.1)),
  yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "",
               range = c(0.5, n_tips + 0.5), autorange = "reversed")
)


# ТЕПЛОВАЯ КАРТА (правая часть)

hover_text <- matrix(
  paste0("<b>Taxon 1:</b> ", rownames(probmat)[row(probmat)],
         "<br><b>Taxon 2:</b> ", colnames(probmat)[col(probmat)],
         "<br><b>PP(conspecific):</b> ", sprintf("%.2f", probmat)),
  nrow = n_tips)

custom_colorscale <- list(
  list(0.00, "#F0FFFF"),  # Почти белый
  list(0.25, "#BDECB6"),  # Светло-зелёный
  list(0.50, "#98FB98"),  # Средне-зелёный
  list(0.95, "#34C924"),  # Тёмно-зелёный
  list(1.00, "#0A5F38")   # Тёмно-тёмно-зелёный (максимальная уверенность)
)

fig_heat <- plotly::plot_ly(
  z = probmat,
  x = 1:n_tips,
  y = 1:n_tips,
  type = "heatmap",
  colorscale = custom_colorscale,  # Или встроенная палитра viridis
  zmin = 0, zmax = 1,
  text = hover_text, hoverinfo = "text",
  showscale = TRUE,
  colorbar = list(
    title = "PP", 
    len = 0.5, 
    x = 1.02,
    tickvals = c(0, 0.25, 0.5, 0.75, 1),
    ticktext = c("0%", "25%", "50%", "75%", "100%")
  )
)

# Настройка осей матрицы (совпадает с деревом)
fig_heat <- fig_heat %>% plotly::layout(
  xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "",
               range = c(0.5, n_tips + 0.5)),
  yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "",
               range = c(0.5, n_tips + 0.5), autorange = "reversed")
)


# ОБЪЕДИНЕНИЕ ЧЕРЕЗ SUBPLOT (1:1)

fig_combined <- plotly::subplot(
  fig_tree, fig_heat,
  nrows = 1,
  widths = c(0.5, 0.5), # Строго 1:1
  shareY = TRUE,         # Жесткая синхронизация оси Y
  titleX = FALSE, titleY = FALSE
)

fig_combined <- fig_combined %>% plotly::layout(
  title = list(text = "bGMYC4 Conspecificity Probabilities + Consensus Tree", x = 0.5),
  plot_bgcolor = "white",
  hovermode = "closest",
  margin = list(l = 10, r = 60, t = 80, b = 10)
)

# Сохранение
html_path <- file.path(DELIM_DIR, "bGMYC_interactive_heatmap.html")
htmlwidgets::saveWidget(fig_combined, html_path, selfcontained = TRUE, title = "bGMYC Interactive Heatmap")
cat(sprintf("✅ Сохранено: %s\n", normalizePath(html_path)))
cat("🌐 Откройте HTML в браузере. Дерево слева, матрица справа (1:1), всё выровнено по Y.\n\n")


# Функция конвертации списка кластеров в датафрейм
create_motu_df <- function(out_list, method_name) {
  motu_list <- lapply(seq_along(out_list), function(i) {
    data.frame(
      Sequence = out_list[[i]],
      MOTU = paste0(method_name, "_", i),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, motu_list)
}

# Автоматическая генерация таблиц
cat("\n💾 Сохраняю таблицы делимитации...\n")

# Стандартный порог (p = 0.05)
out_005 <- bgmyc.point(probmat, ppcutoff = 0.05)
df_005 <- data.frame(
  Sequence = unlist(out_005),
  MOTU_bGMYC = rep(seq_along(out_005), lengths(out_005)),
  stringsAsFactors = FALSE
)
write.table(df_005, file = file.path(DELIM_DIR, "Delimitation_bGMYC_005.csv"), 
            row.names = FALSE, sep = ";", dec = ".", quote = FALSE, fileEncoding = "UTF-8")
cat(sprintf("   ✅ p=0.05: Delimitation_bGMYC_005.csv (%d кластеров)\n", length(out_005)))

# Строгий порог (p = 0.01)
out_001 <- bgmyc.point(probmat, ppcutoff = 0.01)
df_001 <- data.frame(
  Sequence = unlist(out_001),
  MOTU_bGMYC = rep(seq_along(out_001), lengths(out_001)),
  stringsAsFactors = FALSE
)
write.table(df_001, file = file.path(DELIM_DIR, "Delimitation_bGMYC_001.csv"), 
            row.names = FALSE, sep = ";", dec = ".", quote = FALSE, fileEncoding = "UTF-8")
cat(sprintf("   ✅ p=0.01: Delimitation_bGMYC_001.csv (%d кластеров)\n", length(out_001)))

# Экспорт полной таблицы вероятностей
spec_out  <- bgmyc.spec(final_res)
write.csv(spec_out$specprobs, file.path(DELIM_DIR, "bGMYC_delimitation_results.csv"), 
          row.names = FALSE, fileEncoding = "UTF-8")
cat("   📊 Вероятности: bGMYC_delimitation_results.csv\n")

cat("\n📋 Топ-10 кластеров (p=0.05):\n")
for (i in seq_along(out_005)[1:min(10, length(out_005))]) {
  taxa  <- out_005[[i]]
  cat(sprintf("   %2d: %s%s\n", i, paste(taxa[1:min(3, length(taxa))], collapse = ", "),
              if (length(taxa) > 3) "..." else ""))
}

cat("\n🎉 Интерактивный анализ bGMYC4 завершён!\n")