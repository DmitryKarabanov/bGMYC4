#═══════════════════════════════════════════════════════════════
#bGMYC4 v4.3.0: Полный интерактивный пайплайн
#(параметры -> тест -> мульти/сингл -> визуализация с posterior-colored деревом)
#═══════════════════════════════════════════════════════════════

library(mcmcse)
library(plotly)
library(base64enc)
library(htmlwidgets)
library(treeio)
library(dplyr)
library(ggtree)
library(ape)
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

#═══════════════════════════════════════════════════════════════
# 7.5. НАСТОЯЩИЙ Gelman-Rubin R̂ (несколько цепей на ОДНОМ дереве)
#═══════════════════════════════════════════════════════════════

n_chains_run <- 4L

cat(sprintf("\n📐 Запуск %d независимых цепей (dispersed start) для R̂...\n", n_chains_run))

chain_runs <- bgmyc.chains(
  phylo = tree,
  mcmc = params$mcmc,
  burnin = params$burnin,
  thinning = params$thinning,
  py1 = params$py1, py2 = params$py2,
  pc1 = params$pc1, pc2 = params$pc2,
  t1 = params$t1, t2 = params$t2,
  scale = params$scale,
  start = params$start,
  nchains = n_chains_run,
  seed = 42,
  cores = max(1L, min(n_chains_run, parallel::detectCores(logical = FALSE) - 1L)),
  dispersed = TRUE
)

conv <- bgmyc.convergence(chain_runs)

cat("\n=== Краткая сводка bgmyc.convergence() ===\n")
print(conv$summary)

cat("\n=== R̂ по параметрам ===\n")
print(round(conv$rhat, 3))

cat("\n=== ESS по параметрам ===\n")
print(round(conv$ess_total, 0))

cat("\n=== Acceptance rates по цепям ===\n")
print(round(conv$per_chain_accept, 3))

# Проверка целей
rhat_ok <- all(conv$rhat <= 1.05, na.rm = TRUE)
ess_ok  <- all(conv$ess_total >= 200, na.rm = TRUE)

cat("\n=== Итог диагностики ===\n")
if (rhat_ok) {
  cat("✅ Split-R̂: все значения <= 1.05.\n")
} else {
  bad_rhat <- names(conv$rhat)[conv$rhat > 1.05]
  cat(sprintf("⚠️  Split-R̂ > 1.05 у параметров: %s\n",
              paste(bad_rhat, collapse = ", ")))
  cat("   Увеличьте mcmc/burnin или подстройте scale.\n")
}

if (ess_ok) {
  cat("✅ ESS: все значения >= 200.\n")
} else {
  bad_ess <- names(conv$ess_total)[conv$ess_total < 200]
  cat(sprintf("⚠️  ESS < 200 у параметров: %s\n",
              paste(bad_ess, collapse = ", ")))
  cat("   Увеличьте mcmc или уменьшите thinning.\n")
}

cat("\nℹ️  Число цепей смотрите не по length(chain_runs), а по nrow(conv$per_chain_ess).\n")
cat(sprintf("   Цепей фактически: %d\n", nrow(conv$per_chain_ess)))

# 8. ПОДГОТОВКА ФИНАЛЬНЫХ РЕЗУЛЬТАТОВ
if (analysis_mode == 1) {
  cat("\n🌲 Анализ на одном дереве завершён. Формирую выводы...\n")
  final_res <- list(res_single)
  class(final_res) <- "multibgmyc"
} else {
  cat("\n⚡ Запуск bgmyc.multiphylo на выбранных деревьях...\n")
  n_physical <- parallel::detectCores(logical = FALSE)
  n_workers <- max(1L, min(n_physical - 1, n_sample))
  cat(sprintf("   Воркеры: %d (физических ядер: %d)\n", n_workers, n_physical))

  final_res <- bgmyc.multiphylo(
    multiphylo = trees_sample,
    mcmc = params$mcmc, burnin = params$burnin, thinning = params$thinning,
    py1 = params$py1, py2 = params$py2, pc1 = params$pc1, pc2 = params$pc2,
    t1 = params$t1, t2 = params$t2, scale = params$scale, start = params$start,
    cores = n_workers, seed = 42
  )
  cat("✅ Все деревья успешно обработаны.\n")

  # 📐 Per-tree diagnostics / Диагностика по каждому дереву
  cat("\n📐 Диагностика по каждому дереву (per-tree)...\n")
  acc <- t(vapply(final_res, function(r) as.numeric(r$accept)[1:3], numeric(3)))
  colnames(acc) <- c("py", "pc", "th")
  cat("   Acceptance rates (медиана по деревьям): ",
      paste(sprintf("%s=%.3f", colnames(acc), apply(acc, 2, median)), collapse = " | "),
      "\n", sep = "")

  if (requireNamespace("mcmcse", quietly = TRUE)) {
    ess_th <- vapply(final_res, function(r) {
      e <- tryCatch(as.numeric(mcmcse::ess(r$par[, 3])), error = function(e) NA_real_)
      if (length(e) != 1L || !is.finite(e)) NA_real_ else e
    }, numeric(1))
    if (any(is.finite(ess_th))) {
      cat(sprintf("   ESS(threshold): min=%d | median=%d | max=%d  [цель: >200]\n",
                  round(min(ess_th, na.rm = TRUE)),
                  round(median(ess_th, na.rm = TRUE)),
                  round(max(ess_th, na.rm = TRUE))))
      bad <- which(is.finite(ess_th) & ess_th < 200)
      if (length(bad)) {
        cat(sprintf("   ⚠️  ESS < 200 у %d дерев(ьев): %s — увеличьте mcmc или уменьшите thinning.\n",
                    length(bad), paste(head(bad, 10), collapse = ", ")))
      } else {
        cat("   ✅ ESS достаточен у всех деревьев.\n")
      }
    } else {
      cat("   ⚠️  ESS рассчитать не удалось (цепи слишком короткие).\n")
    }
  } else {
    cat("   ⚠️  ESS: пакет 'mcmcse' не установлен — install.packages('mcmcse')\n")
  }

}


#═══════════════════════════════════════════════════════════════
# 9. ВИЗУАЛИЗАЦИЯ И АВТОМАТИЧЕСКИЙ ЭКСПОРТ (КАК В РАБОЧЕМ ПРИМЕРЕ)
#═══════════════════════════════════════════════════════════════
cat("\n🌡️ Построение интерактивной карты вероятностей...\n")
probmat <- spec.probmat(final_res)

# ── Округление матрицы ──
probmat <- round(probmat, 4)

p_tree <- suppressWarnings(ggtree(tree, layout = "rectangular"))
tips_data <- p_tree$data %>% dplyr::filter(isTip) %>% dplyr::arrange(y)
tip_order <- tips_data$label
n_tips <- length(tip_order)

if (!all(rownames(probmat) == tip_order) || !all(colnames(probmat) == tip_order)) {
  probmat <- probmat[tip_order, tip_order]
}

cat("\n🎨 Создание кастомной визуализации с деревом (1:1)...\n")

if (!exists("tree_beast") || !inherits(tree_beast, "treedata")) {
  tree_beast <- treeio::read.beast(consensus_path)
}

if(!is.null(tree$edge.length)) tree$edge.length[tree$edge.length < 0] <- 0
options(ignore.negative.edge = TRUE)
p_tree_raw <- suppressWarnings(ggtree(tree, layout = "rectangular"))
tree_data_raw <- p_tree_raw$data
tips_data_raw <- tree_data_raw %>% dplyr::filter(isTip) %>% dplyr::arrange(y)
tip_order <- tips_data_raw$label

common_taxa <- intersect(tip_order, rownames(probmat))
if (length(common_taxa) == 0) stop("❌ Нет общих таксонов между деревом и матрицей!")

tips_to_drop <- setdiff(tree_beast@phylo$tip.label, common_taxa)
if(length(tips_to_drop) > 0) {
  tree_beast_filtered <- treeio::drop.tip(tree_beast, tips_to_drop)
} else {
  tree_beast_filtered <- tree_beast
}

p_tree <- suppressWarnings(ggtree(tree_beast_filtered, layout = "rectangular"))
tree_data <- p_tree$data
tips_data <- tree_data %>% dplyr::filter(isTip) %>% dplyr::arrange(y)
tip_order <- tips_data$label
n_tips <- length(tip_order)

probmat <- probmat[tip_order, tip_order]
cat(sprintf("   ✅ Матрица %d×%d синхронизирована с деревом\n", n_tips, n_tips))

# Извлекаем posterior
post_col <- intersect(c("posterior", "prob", "Posterior", "PROB"), colnames(tree_data))
if(length(post_col) > 0) {
  post_col <- post_col[1]
  node_posterior <- tree_data[[post_col]]
  names(node_posterior) <- tree_data$node
} else {
  node_posterior <- sapply(tree_data$label, function(lbl) {
    if(is.na(lbl) || lbl == " ") return(NA_real_)
    m <- regmatches(lbl, regexpr("(posterior|prob)\\s*=\\s*([0-9\\.eE\\-]+)", lbl, ignore.case = TRUE, perl = TRUE))
    if(length(m) == 0) return(NA_real_)
    num_str <- sub("^(posterior|prob)\\s*=\\s*", "", m, ignore.case = TRUE, perl = TRUE)
    val <- as.numeric(num_str)
    if(!is.na(val) && val > 1) val <- val / 100
    return(val)
  })
  names(node_posterior) <- tree_data$node
}

edges <- tree_data %>% dplyr::filter(!is.na(parent))

get_pp_color <- function(pp) {
  if(is.na(pp)) return("#CCCCCC")
  pp <- max(0, min(1, pp))
  colors <- list(
    c(0.0, 1.0, 0.0, 0.0), c(0.25, 1.0, 0.5, 0.0), c(0.5, 1.0, 1.0, 0.0),
    c(0.75, 0.5, 1.0, 0.0), c(1.0, 0.0, 0.7, 0.0)
  )
  if(pp <= 0) return(sprintf("#%02X%02X%02X", 255, 0, 0))
  if(pp >= 1) return(sprintf("#%02X%02X%02X", 0, 179, 0))
  for(i in 1:(length(colors)-1)) {
    if(pp >= colors[[i]][1] && pp <= colors[[i+1]][1]) {
      t <- (pp - colors[[i]][1]) / (colors[[i+1]][1] - colors[[i]][1])
      r <- colors[[i]][2] + t * (colors[[i+1]][2] - colors[[i]][2])
      g <- colors[[i]][3] + t * (colors[[i+1]][3] - colors[[i]][3])
      b <- colors[[i]][4] + t * (colors[[i+1]][4] - colors[[i]][4])
      return(sprintf("#%02X%02X%02X", round(r*255), round(g*255), round(b*255)))
    }
  }
  return("#CCCCCC")
}

# Цветовая шкала для ТЕПЛОВОЙ КАРТЫ (PP conspecific) — бело-зелёная
custom_colorscale_heat <- list(
  c(0.00, "#F0FFFF"),  # Почти белый (низкая вероятность)
  c(0.25, "#BDECB6"),  # Светло-зелёный
  c(0.50, "#98FB98"),  # Средне-зелёный
  c(0.95, "#34C924"),  # Тёмно-зелёный
  c(1.00, "#0A5F38")   # Тёмно-тёмно-зелёный (максимальная уверенность)
)

# ── ДЕРЕВО (левая часть) ──
parent_x_lookup <- setNames(tree_data$x, tree_data$node)
parent_y_lookup <- setNames(tree_data$y, tree_data$node)
px <- parent_x_lookup[as.character(edges$parent)]
py <- parent_y_lookup[as.character(edges$parent)]
cx <- round(edges$x, 4); cy <- round(edges$y, 4)
px <- round(px, 4); py <- round(py, 4)

child_nodes <- as.character(edges$node)
pp_vals <- node_posterior[child_nodes]
branch_colors <- sapply(pp_vals, get_pp_color)
pp_text <- ifelse(is.na(pp_vals), "N/A", sprintf("%.2f", pp_vals))
hover_txt <- paste0("<b>Node: </b>", edges$node, "<br><b>Posterior: </b>", pp_text)

fig_tree <- plotly::plot_ly() %>%
  plotly::add_segments(
    data = data.frame(x=px, xend=cx, y=cy, yend=cy, col=branch_colors, hov=hover_txt, stringsAsFactors = FALSE),
    x=~x, xend=~xend, y=~y, yend=~yend, line=list(width=5, color=~I(col)),
    hovertext=~hov, hoverinfo="text", showlegend=FALSE
  ) %>%
  plotly::add_segments(
    data = data.frame(x=px, xend=px, y=py, yend=cy, col=branch_colors, hov=hover_txt, stringsAsFactors = FALSE),
    x=~x, xend=~xend, y=~y, yend=~yend, line=list(width=5, color=~I(col)),
    hovertext=~hov, hoverinfo="text", showlegend=FALSE
  )

max_x <- max(tree_data$x, na.rm = TRUE)
label_offset <- max_x * 0.05
label_step <- max(1, ceiling(n_tips / 40))
tips_labeled <- tips_data %>% dplyr::mutate(row_num = dplyr::row_number()) %>% dplyr::filter(row_num %% label_step == 1)

fig_tree <- fig_tree %>% plotly::add_text(
  data = tips_labeled, x = ~round(x + label_offset, 4), y = ~round(y, 4), text = ~label,
  textposition = "middle left", textfont = list(size = 16, family = "monospace", color = "#333333"),
  showlegend = FALSE, hoverinfo = "skip") %>% plotly::layout(
    xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "",
                 range = c(-0.02, round(max_x + label_offset + max_x*0.1, 4))),
    yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "",
                 range = c(0.5, n_tips + 0.5), autorange = "reversed")
  )

# ── ТЕПЛОВАЯ КАРТА (КАТЕГОРИАЛЬНЫЕ ОСИ КАК В РАБОЧЕМ ПРИМЕРЕ!) ──
fig_heat <- plotly::plot_ly(
  z = probmat,
  x = tip_order,  # ← Категориальная ось!
  y = tip_order,  # ← Категориальная ось!
  type = "heatmap",
  colorscale = custom_colorscale_heat,
  zmin = 0, zmax = 1,
  hovertemplate = paste0(
    "<b>Taxon 1: </b>%{y}<br>",
    "<b>Taxon 2: </b>%{x}<br>",
    "<b>PP(conspecific): </b>%{z:.2f}",
    "<extra></extra>"
  ),
  showscale = TRUE,
  colorbar = list(
    title = "PP",
    len = 0.5,
    x = 1.02,
    tickvals = c(0, 0.25, 0.5, 0.75, 1),
    ticktext = c("0%", "25%", "50%", "75%", "100%")
  )
) %>% plotly::layout(
  xaxis = list(title = "", showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
  yaxis = list(title = "", showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE, autorange = "reversed")
)

# ── ОБЪЕДИНЕНИЕ (КАК В РАБОЧЕМ ПРИМЕРЕ!) ──
fig_combined <- plotly::subplot(
  fig_tree, fig_heat,
  nrows = 1,
  widths = c(0.50, 0.50),
  shareY = FALSE,  # ← ВАЖНО: убираем, оси разного типа
  titleX = FALSE, titleY = FALSE
) %>% plotly::layout(
  yaxis = list(domain = c(0, 1)),   # Ось дерева
  yaxis2 = list(domain = c(0, 1)),  # Ось матрицы
  title = list(text = "bGMYC4 Conspecificity Probabilities + Consensus Tree", x = 0.5),
  plot_bgcolor = "white",
  margin = list(l = 10, r = 80, t = 80, b = 30)
)

# ── ЭКСПОРТ HTML ──
html_path <- file.path(DELIM_DIR, "bGMYC_interactive_heatmap.html")
htmlwidgets::saveWidget(fig_combined, html_path, selfcontained = FALSE, title = "bGMYC Interactive Heatmap")

cat(sprintf("✅ Сохранено: %s\n", normalizePath(html_path)))
cat(sprintf("📦 Размер файла: %.1f MB\n", file.info(html_path)$size / 1024^2))
cat("🌐 Откройте HTML. Hover показывает имена таксонов + PP (стандартный Plotly)!\n\n")

# # EN: Export clusters at different thresholds / Экспорт кластеров на разных порогах
thresholds <- c(0.95, 0.5, 0.0, 0.01, 0.001)
for (p in thresholds) {
  out <- bgmyc.point(probmat, ppcutoff = p)
  df  <- data.frame(
    Sequence   = unlist(out),
    MOTU_bGMYC = rep(seq_along(out), lengths(out)),
    stringsAsFactors = FALSE
  )
  fname <- file.path(DELIM_DIR, sprintf("Delimitation_bGMYC_%.3f.csv", p))
  utils::write.table(df, file = fname, row.names = FALSE,
                     sep = ";", dec = ".", quote = FALSE)
  cat(sprintf("✅ ppcutoff=%.3f → %d кластеров, %d строк → %s\n",
              p, length(out), nrow(df), basename(fname)))
}

spec_out <- bgmyc.spec(final_res)
utils::write.csv(spec_out$specprobs,
                 file = file.path(DELIM_DIR, "bGMYC_full_probs.csv"),
                 row.names = FALSE)
cat("✅ bGMYC_full_probs.csv сохранён в", DELIM_DIR, "\n")

cat("\n🎉 Интерактивный анализ bGMYC4 завершён!\n")
