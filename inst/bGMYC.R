# ═══════════════════════════════════════════════════════════════
# bGMYC4 v4.1.0: Полный интерактивный пайплайн (параметры -> тест -> мульти/сингл)
# ═══════════════════════════════════════════════════════════════
library(ape)
library(future)
library(future.apply)
library(mcmcse)
library(plotly)
library(base64enc)
library(htmlwidgets)
library(bGMYC4)

cat("✅ Все критичные зависимости загружены\n")

cat("=== 🌿 bGMYC4 Interactive Analysis ===\n\n")

# ── 1. ЗАГРУЗКА ДЕРЕВЬЕВ ─────────────────────────────────────
consensus_path <- readline(prompt = "📥 Путь к консенсусному дереву (tree): ")
posterior_path <- readline(prompt = "📥 Путь к набору деревьев (trees): ")
if (!file.exists(consensus_path)) stop("❌ Консенсусное дерево не найдено!")
if (!file.exists(posterior_path)) stop("❌ Файл с деревьями не найден!")

tree <- read.nexus(consensus_path)
all_trees <- read.nexus(posterior_path)
class(all_trees) <- "multiPhylo"
cat(sprintf("✅ Загружено: 1 консенсус + %d деревьев в апостериоре\n", length(all_trees)))

# ── 2. ВЫБОР РЕЖИМА АНАЛИЗА ──────────────────────────────────
cat("\n📊 Выберите режим анализа:\n")
cat("   1. Только консенсусное дерево (быстрый анализ, без учёта топологической неопределённости)\n")
cat("   2. Консенсус + множество деревьев (учёт филогенетической неопределённости)\n")
mode_str <- readline("   Введите 1 или 2 (по умолчанию 2): ")
analysis_mode <- ifelse(nchar(trimws(mode_str)) == 0 || is.na(as.integer(mode_str)), 2, as.integer(mode_str))
if (!analysis_mode %in% c(1, 2)) {
  cat("⚠️ Неверный выбор. Установлен режим 2.\n")
  analysis_mode <- 2
}

# ── 3. ВЫБОРКА ДЕРЕВЬЕВ (ТОЛЬКО ДЛЯ РЕЖИМА 2) ────────────────
if (analysis_mode == 2) {
  n_total <- length(all_trees)
  burnin_idx <- floor(n_total * 0.10)
  n_str <- readline(prompt = sprintf("\n🎲 Сколько деревьев выбрать? (по умолчанию 10, доступно: %d): ", n_total - burnin_idx))
  n_sample <- ifelse(nchar(n_str) == 0 || is.na(as.numeric(n_str)), 10, as.integer(n_str))
  if (n_sample > (n_total - burnin_idx)) n_sample <- n_total - burnin_idx
  
  set.seed(42) # Воспроизводимость
  trees_sample <- all_trees[sample((burnin_idx + 1):n_total, n_sample)]
  class(trees_sample) <- "multiPhylo"
  cat(sprintf("✅ Выбрано %d случайных деревьев (индексы %d–%d)\n", n_sample, burnin_idx + 1, n_total))
} else {
  cat("\n📦 Режим 1 выбран: анализ только на консенсусном дереве.\n")
  trees_sample <- NULL # Не занимаем память
}

# ── 4. ПРОВЕРКА УЛЬТРАМЕТРИЧНОСТИ ────────────────────────────
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

# ── 5. АУТГРУППА (ПОДДЕРЖКА НЕСКОЛЬКИХ ТАКСОНОВ) ─────────────
og_prompt <- readline(prompt = "\n🌳 Удалять аутгруппы? (y/n, по умолчанию n): ")
if (tolower(og_prompt) == "y") {
  og_str <- readline(prompt = "   Введите имена таксонов-аутгрупп (через запятую): ")
  og_names <- trimws(unlist(strsplit(og_str, "[,;\\s]+")))
  og_names <- og_names[nchar(og_names) > 0]
  
  if (length(og_names) > 0) {
    existing <- og_names[og_names %in% tree$tip.label]
    missing  <- setdiff(og_names, existing)
    
    if (length(missing) > 0) {
      cat(sprintf("⚠️  Не найдены в дереве (будут пропущены): %s\n", paste(missing, collapse = ", ")))
    }
    if (length(existing) > 0) {
      tree <- drop.tip(tree, existing)
      if (analysis_mode == 2) {
        trees_sample <- lapply(trees_sample, drop.tip, existing)
        class(trees_sample) <- "multiPhylo"
      }
      cat(sprintf("✅ Удалено %d таксонов: %s\n", length(existing), paste(existing, collapse = ", ")))
    } else {
      cat("⚠️  Ни один из указанных таксонов не найден. Пропускаю.\n")
    }
  }
}
ntips <- length(tree$tip.label)
cat(sprintf("📊 Итоговый размер: %d таксонов\n\n", ntips))

# ── 6. БЛОК ЗАДАНИЯ ПАРАМЕТРОВ (ВЫПОЛНЯЕТСЯ ДО ПЕРВОГО ЗАПУСКА) ─
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
params$mcmc     <- input_scalar("mcmc", 10000)
params$burnin   <- input_scalar("burnin", 1000)
params$thinning <- input_scalar("thinning", 10)
params$py1      <- input_scalar("py1", 0)
params$py2      <- input_scalar("py2", 1.5)
params$pc1      <- input_scalar("pc1", 0)
params$pc2      <- input_scalar("pc2", 2)
params$t1       <- input_scalar("t1", 2)
params$t2       <- input_scalar("t2", min(ntips - 1, 100))
params$scale    <- input_vector("scale", c(20, 10, 5))
params$start    <- input_vector("start", c(1, 1, floor((params$t1 + min(ntips-1, 50))/2)))

# Валидация границ априоров и старта
params$t2 <- max(params$t1 + 2, min(params$t2, ntips - 1))
params$start[3] <- max(params$t1 + 1, min(params$start[3], params$t2 - 1))
cat(sprintf("🔒 Параметры зафиксированы: t ∈ [%d, %d] | start[3] = %d\n\n", params$t1, params$t2, params$start[3]))

# ── 7. ИНТЕРАКТИВНЫЙ ЦИКЛ ТЕСТА SINGLEPHY ────────────────────
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
  
  # ── 🟢 ДИНАМИЧЕСКАЯ АНАЛИТИКА СХОДИМОСТИ ────────────────────
  ar <- res_single$accept
  cat(sprintf("\n📊 Acceptance rates: py=%.3f | pc=%.3f | th=%.3f\n", ar[1], ar[2], ar[3]))
  
  # 🔢 ESS Calculation (Effective Sample Size)
  if (requireNamespace("mcmcse", quietly = TRUE)) {
    ess_vals <- sapply(1:4, function(col) round(mcmcse::ess(res_single$par[, col])))
    cat(sprintf("🔢 ESS: py=%d | pc=%d | th=%d | logL=%d [Оптимум: >200]\n",
                ess_vals[1], ess_vals[2], ess_vals[3], ess_vals[4]))
  } else {
    cat("⚠️  ESS: пакет 'mcmcse' не установлен. Выполните install.packages('mcmcse')\n")
    ess_vals <- NULL
  }
  
  recs  <- character(0)
  param_names  <- c("py", "pc", "th")
  for (i in 1:3) {
    if (ar[i]  > 0.55) {
      new_val  <- round(params$scale[i] * 1.5)
      if (new_val == params$scale[i]) new_val  <- params$scale[i] + (if(i  < 3) 5 else 2)
      recs  <- c(recs, sprintf("• %s: слишком высокая (%.2f) → увеличьте scale[%d] с %g до ~%d",
                               param_names[i], ar[i], i, params$scale[i], new_val))
    } else if (ar[i]  < 0.15) {
      new_val  <- max(if(i  < 3) 2 else 1, round(params$scale[i] * 0.5))
      recs  <- c(recs, sprintf("• %s: слишком низкая (%.2f) → уменьшите scale[%d] с %g до ~%d",
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
  # ─────────────────────────────────────────────────────────────
  
  next_prompt <- if (analysis_mode == 1) {
    "⏭️  Перейти к финальному выводу (y) или изменить параметры (n)? [y/n]: "
  } else {
    "⏭️  Перейти к multiphylo (y) или изменить параметры (n)? [y/n]: "
  }
  
  choice <- readline(prompt = sprintf("\n%s", next_prompt))
  if (tolower(choice) != "n") break
  
  cat("\n📝 Обновите параметры (Enter = оставить текущее):\n")
  params$mcmc     <- input_scalar("mcmc", params$mcmc)
  params$burnin   <- input_scalar("burnin", params$burnin)
  params$thinning <- input_scalar("thinning", params$thinning)
  params$py1      <- input_scalar("py1", params$py1)
  params$py2      <- input_scalar("py2", params$py2)
  params$pc1      <- input_scalar("pc1", params$pc1)
  params$pc2      <- input_scalar("pc2", params$pc2)
  params$t1       <- input_scalar("t1", params$t1)
  params$t2       <- input_scalar("t2", params$t2)
  params$scale    <- input_vector("scale", params$scale)
  params$start    <- input_vector("start", params$start)
  
  if (length(params$scale) != 3) params$scale <- c(20, 10, 5)
  if (length(params$start) != 3) params$start <- c(1, 0.5, floor((params$t1 + params$t2)/2))
  params$t2 <- max(params$t1 + 2, min(params$t2, ntips - 1))
  params$start[3] <- max(params$t1 + 1, min(params$start[3], params$t2 - 1))
  cat(sprintf("🔒 Параметры обновлены: t ∈ [%d, %d] | start[3] = %d\n\n", params$t1, params$t2, params$start[3]))
}

# ── 8. ПОДГОТОВКА ФИНАЛЬНЫХ РЕЗУЛЬТАТОВ (УНИВЕРСАЛЬНЫЙ ОБЪЕКТ) ─
if (analysis_mode == 1) {
  cat("\n🌲 Анализ на одном дереве завершён. Формирую выводы...\n")
  final_res <- list(res_single)
  class(final_res) <- "multibgmyc"
} else {
  cat("\n⚡ Запуск bgmyc.multiphylo на выбранных деревьях...\n")
  n_physical <- parallel::detectCores(logical = FALSE)
  n_workers  <- min(n_physical - 1, n_sample)
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
  
  # 📐 Gelman-Rubin (R̂) диагностика по цепям разных деревьев
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

# ─── 9. ИНТЕРАКТИВНАЯ ТЕПЛОВАЯ КАРТА И ДЕЛИМИТАЦИЯ ──────────────────────────
cat("\n🌡️ Построение интерактивной матрицы конспецифичности...\n")
probmat <- spec.probmat(final_res)

# 🎨 Запрос параметров визуализации
pal_prompt <- readline(prompt = "🎨 Выберите палитру (green/viridis/RdYlBu/classic) [по умолчанию green]: ")
palette <- ifelse(nchar(trimws(pal_prompt)) == 0 || !trimws(pal_prompt) %in% c("green", "viridis", "RdYlBu", "classic"), "green", trimws(pal_prompt))

tw_prompt <- readline(prompt = "🌲 Ширина дерева (0.1–0.5) [по умолчанию 0.25]: ")
tree_width <- suppressWarnings(as.numeric(tw_prompt))
if (is.na(tree_width) || tree_width < 0.1 || tree_width > 0.5) tree_width <- 0.25

# 🌐 Генерация HTML
html_path <- "bGMYC_interactive_heatmap.html"
plot.interactive.probmat(
  x = probmat, tree = tree, palette = palette, tree_width = tree_width, save_html = html_path
)
cat(sprintf("✅ Интерактивный график сохранён: %s (палитра: %s)\n", normalizePath(html_path), palette))
cat("🌐 Откройте файл в браузере для зума, перемещения и просмотра вероятностей.\n\n")

# Продолжаем без ожидания
ppc_str <- readline(prompt = "🔪 Введите ppcutoff для выделения видов (по умолчанию 0.05): ")
ppcutoff <- ifelse(nchar(ppc_str) == 0 || is.na(as.numeric(ppc_str)), 0.05, as.numeric(ppc_str))
out <- bgmyc.point(probmat, ppcutoff = ppcutoff)
cat(sprintf("\n🔍 Выделено кластеров при PP > %.2f: %d\n", ppcutoff, length(out)))

# Экспорт
spec_out <- bgmyc.spec(final_res)
write.csv(spec_out$specprobs, "bGMYC_delimitation_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("📊 Таблица вероятностей сохранена: bGMYC_delimitation_results.csv\n")
cat("\n📋 Топ-10 кластеров по вероятности:\n")
for (i in seq_along(out)[1:min(10, length(out))]) {
  taxa <- out[[i]]
  cat(sprintf("   %2d: %s%s\n", i, paste(taxa[1:min(3, length(taxa))], collapse = ", "),
              if (length(taxa) > 3) "..." else ""))
}

cat("📋Полный результат делимитации:")
print(out)

# Превращаем список 'out' в таблицу и сохраняем
df <- data.frame(
  Sequence   = unlist(out),
  MOTU_bGMYC = rep(seq_along(out), sapply(out, length))
)

write.table(df, 
            file = "C:/GEN/Bosmina/Delimitation/Delimitation_bGMYC.csv",
            sep = ";", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")

cat(sprintf("✅ Сохранено: %d таксонов, %d MOTU\n", nrow(df), length(out)))

cat("\n🎉 Интерактивный анализ bGMYC4 завершён!\n")
