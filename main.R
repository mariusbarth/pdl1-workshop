## ----setup, include = FALSE----------------------------------------------
# load packages
library("papaja")
library("afex")
library("emmeans")


## ----analysis-preferences------------------------------------------------
# Seed for random number generation
set.seed(42)
knitr::opts_chunk$set(
  cache.extra = knitr::rand_seed
  , message = FALSE
  , warning = FALSE
)


## ----prepare-data--------------------------------------------------------
acquisition <- readRDS("data/acquisition-task.rds")
generation <- readRDS("data/generation-task.rds")

# set some variable labels for pretty plotting
variable_labels(acquisition) <- c(
  age = "Age"
  , sex = "Sex"
  , material = "Material"
  , generation = "Generation task"
  , order = "Block order"
  , block = "Block number"
  , trial = "Trial number"
  , S = "Stimulus location"
  , R = "Response location"
  , SRI = "Response time [ms]"
  , error = "Error"
  , frequency = "Location frequency"
)

variable_labels(generation) <- c(
  age = "Age"
  , sex = "Sex"
  , material = "Material"
  , generation = "Generation task"
  , order = "Block order"
  , block = "Block number"
  , trial = "Trial number"
  , PD_instruction = "PD instruction"
  , correct_SOC = "Correctly generated second-order conditional"
)


## ----participants--------------------------------------------------------
# stats for participants section
excluded_participants <- unique(c(55, 74, 75, 126, 5, 7, 9, 12, 15))

N <- length(unique(generation$id))

n.excluded_participants <- length(excluded_participants)

tmp <- aggregate(formula = correct_SOC ~ id + sex + age, data = generation, FUN=mean)
Sex <- table(tmp[["sex"]])

meanAge <- printnum(mean(tmp$age), digits = 1)
rangeAge <- paste(range(tmp$age), collapse = "-")


## ----acquisition-rt, fig.cap = "(ref:acquisition-rt)"--------------------
# Reaction times during SRTT phase (i.e., training or acquisition)
tmp <- acquisition[
  acquisition$error == 0 &
  acquisition$material != "No-learning" & # these Ss didn't work on SRTT
  acquisition$trial > 1 &
  acquisition$included_participant,
]

out <- apa_print(
  aov_ez(data = tmp, dv = "SRI", id = "id", between = "material", within = "block")
)

par(mfrow = c(1, 2)) # place two different plots side-by-side

# left panel: RTs in permuted and random material groups
apa_lineplot(
  data = tmp
  , id = "id"
  , dv = "SRI"
  , factors = c("block", "material")
  , dispersion = wsci
  , ylab = "Mean RT [ms]"
  , args_lines = list(lty = c("solid", "dotted"))
  , args_points = list(pch = c(21, 23), bg = c("grey70", "white"))
  , args_legend = list(legend = c("Permuted", "Random"))
  , ylim = c(475, 650)
  , jit = .01
)

# within the permuted-material group, high vs. low frequency locations can be distinguished
tmp.perm <- acquisition[
  acquisition$error == 0 &
  acquisition$trial > 1 &
  acquisition$included_participant & 
  acquisition$material == "Permuted" & 
  !is.na(acquisition$frequency),
]

out.perm <- apa_print(
  aov_ez(data = tmp.perm, dv = "SRI", within = c("block", "frequency"), id = "id")
)

apa_lineplot(
  data = tmp.perm
  , id = "id"
  , dv = "SRI"
  , factors = c("block", "frequency")
  , dispersion = wsci
  , ylab = ""
  , args_lines = list(lty = c("solid", "solid"))
  , args_points = list(pch = c(21, 21))
  , args_legend = list(legend = c("High", "Low"), title = "Location frequency")
  , ylim = c(475, 650)
  , jit = .01
)


## ----generation-soc, fig.cap = "(ref:generation-soc)"--------------------
# proportin of correctly generated triplets in the generation task
# excluding repetitions, and trials after repetitions
tmp <- generation[
  generation$included_participant &
  generation$trial > 2 &
  generation$repetition == 0 &
  generation$post_repetition == 0, 
]

# visualize
apa_barplot(
  data = tmp
  , id = "id"
  , dv = "correct_SOC"
  , ylab = "Proportion correct generation"
  , factors = c("material","PD_instruction", "generation")
  , ylim = c(0, 1)
  , main = c("Free generation", "Cued generation")
  , intercept = .2
)

# aggregate once, because many, many ANOVAs will follow
tmp2 <- aggregate(
  formula = correct_SOC ~ id + material + generation + order + PD_instruction
  , data = tmp
  , FUN = mean
)

# ANOVA for complete design
out <- apa_print(
  aov_ez(data = tmp2, id = "id", dv = "correct_SOC", between = c("material", "generation", "order"), within = "PD_instruction")
)

# separate analysis for free-generation task
tmp.a <- tmp2[tmp2$generation == "Free", ]
out.a <- apa_print(
  aov_ez(data = tmp.a, id = "id", dv = "correct_SOC", between = c("material", "order"), within = "PD_instruction")
)

# Post-hoc paired comparisons for free generation group
fit <- aov_ez(data = tmp.a, id = "id", dv = "correct_SOC", between = c("material", "order"), within = "PD_instruction")
tukey <- apa_print(pairs(emmeans(fit, specs = "material"), adjust = "tukey"))



# Analyses for ordinal-PD approach ----
# (you may want to skip these)

# free generation, only exclusion task
tmp.ae <- tmp2[tmp2$generation == "Free" & tmp2$PD_instruction == "Exclusion", ]
out.ae <- apa_print(aov_ez(data = tmp.ae,id = "id", dv = "correct_SOC", between = c("material", "order")))

# free generation, only exclusion task, permuted vs. no-learning condition
tmp.aep0 <- tmp2[tmp2$generation == "Free" & tmp2$PD_instruction == "Exclusion" & tmp2$material!="Random", ]
out.aep0 <- apa_print(aov_ez(data = tmp.aep0, id = "id", dv = "correct_SOC", between = c("material", "order")))

# free generation, only exclusion task, permuted vs. random material
tmp.aepr <- tmp2[tmp2$generation == "Free" & tmp2$PD_instruction == "Exclusion" & tmp2$material!="No-learning", ]
out.aepr <- apa_print(aov_ez(data = tmp.aepr, id = "id", dv = "correct_SOC", between = c("material","order")))

# free generation, exclusion, random vs no-learning condition
tmp.aer0 <- tmp2[tmp2$generation == "Free" & tmp2$PD_instruction == "Exclusion" & tmp2$material!="Permuted", ]
out.aer0 <- apa_print(aov_ez(data = tmp.aer0, id = "id", dv = "correct_SOC", between = c("material","order")))

# free generation, inclusion
tmp.ai<-tmp2[tmp2$generation == "Free" & tmp2$PD_instruction == "Inclusion",]
out.ai<-apa_print(aov_ez(data = tmp.ai, id = "id", dv = "correct_SOC", between = c("material","order")))

# free generation, inclusion, permuted vs. no-learning group
tmp.aip0 <- tmp2[tmp2$generation == "Free" & tmp2$PD_instruction == "Inclusion" & tmp2$material!="Random", ]
out.aip0 <- apa_print(aov_ez(data = tmp.aip0, id = "id", dv = "correct_SOC", between = c("material","order")))

# free generation, inclusion, permuted vs. random group
tmp.aipr <- tmp2[tmp2$generation == "Free" & tmp2$PD_instruction == "Inclusion" & tmp2$material!="No-learning", ]
out.aipr <- apa_print(aov_ez(data = tmp.aipr, id = "id", dv = "correct_SOC", between=c("material","order")))

# free generation, only exclusion task, random vs. no-learning group
tmp.air0 <- tmp2[tmp2$generation == "Free" & tmp2$PD_instruction == "Inclusion" & tmp2$material!="Permuted", ]
out.air0 <- apa_print(aov_ez(data = tmp.air0, id = "id", dv = "correct_SOC", between = c("material","order")))


## ----was-cue-including-reversals, fig.cap = "(ref:was-cue-including-reversals)"----
# Analyze whether participants responded with locations that were also presented
# as a cue on the respective trial (in the cued-generation task)

tmp <- generation[
  generation$included_participant &
  generation[["repetition"]] == 0 &
  generation$generation == "Cued",
]

apa_barplot(
   data = tmp
  , id = "id"
  , dv = "was_cue"
  , ylab = ""
  , factors = c("material","PD_instruction", "block")
  , ylim = c(0, 1)
  , intercept = .6
  , main = c("First block", "Second block")
  , args_legend = list(plot = c(TRUE, FALSE)) # plot legend in left (not right) plot
)

# ANOVA for full design
model <- aov_ez(data = tmp, id = "id", dv = "was_cue", between = c("material", "order"), within = "PD_instruction")
out <- apa_print(model)

# post-hoc pairwise comparisons
tukey <- apa_print(pairs(emmeans(model, specs = "material")))


# analyse blocks separately (extract only `full_result` to output object)
out_1st <- apa_print(
  aov_ez(data = tmp[tmp$block == "first", ], id = "id", dv = "was_cue", between = c("material", "PD_instruction"))
)$full_result

out_2nd <- apa_print(
  aov_ez(data = tmp[tmp$block == "second",], id = "id", dv = "was_cue", between = c("material","PD_instruction"))
)$full_result


## ----was-cue-excluding-reversals, fig.cap = "(ref:was-cue-excluding-reversals)"----
# do the same analyses, but this time exclude reversals
tmp <- generation[
  generation$included_participant &
  generation[["trial"]] > 1 &
  generation[["repetition"]] == 0 &
  generation[["reversal"]] == 0 &
  generation$generation == "Cued"
  ,
]

# ANOVA for full design
out <- apa_print(
  aov_ez(data = tmp, id = "id", dv = "was_cue", between = c("material", "order"), within = "PD_instruction")
)

# separate analyses for first vs. second generation block
out_1st <- apa_print(
  aov_ez(data = tmp[tmp$block=="first", ], id = "id", dv = "was_cue", between = c("material", "PD_instruction"))
)

out_2nd <- apa_print(
  aov_ez(data = tmp[tmp$block=="second", ], id = "id", dv = "was_cue", between = c("material", "PD_instruction"))
)


## ----create_r-references-------------------------------------------------
r_refs(file = "r-references.bib")


## ----echo = FALSE, results = 'asis'--------------------------------------
render_appendix('appendix.rmd')

