# Tests for the shipped analysis skills and their injection into prompts.
#
# The skills exist because the 2026-10-16 sweep measured knowledge failures
# the engine could not prevent: ordering a consumer before its evidence
# producers, passing an OrgDb where the ontology database belongs. These
# tests pin the plumbing (files parse, selection matches, injection reaches
# the prompt); the CONTENT of the skills is documentation, not code, and is
# exercised against a real provider by the benchmark repo, not here.

test_that("skills are discoverable and parse with name, description and body", {
    skills <- AISkills()
    expect_true(length(skills) >= 3L)
    expect_true("sclet-analysis" %in% skills)
    for (s in skills) {
        parsed <- ReadAISkill(s)
        expect_equal(parsed$name, s)
        # a description that failed to parse is the failure that made the
        # first version of the parser useless: descriptions are the
        # selection key
        expect_true(nchar(parsed$description) > 100L,
            info = paste("description too short for", s))
        expect_true(nchar(parsed$body) > 500L, info = s)
        # folded-scalar artifacts must not leak into the text the model sees
        expect_false(grepl("^\\s*>\\s*$", parsed$description))
        expect_false(grepl("\n  ", substr(parsed$description, 1, 200)))
    }
})

test_that("sclet_ai_skill_context selects by topic and formats", {
    all_bodies <- sclet_ai_skill_context()
    expect_true(nchar(all_bodies) > 15000L)
    expect_true(grepl("## Skill: sclet-analysis", all_bodies, fixed = TRUE))

    index <- sclet_ai_skill_context(format = "index")
    expect_false(grepl("## Skill:", index, fixed = TRUE))
    expect_true(grepl("sclet-provenance", index, fixed = TRUE))

    # selection narrows: a claim-only topic must not carry the full
    # analysis body
    only_claims <- sclet_ai_skill_context(c("claim"))
    expect_true(nchar(only_claims) < nchar(all_bodies))
    expect_false(grepl("## Skill: sclet-analysis", only_claims, fixed = TRUE))

    expect_equal(sclet_ai_skill_context(c("zzz_no_match")), "")
})

test_that("goal-to-topic mapping over-includes rather than misses", {
    # over-inclusion is the safe direction: a skill that was not needed
    # costs characters, a missed one costs a repair round
    expect_true("analysis guidelines" %in%
        sclet_ai_goal_skill_topics("Cluster the cells"))
    expect_true("claim" %in%
        sclet_ai_goal_skill_topics("Find markers for each cluster"))
    expect_true("provenance" %in%
        sclet_ai_goal_skill_topics("compare the two analyses and report reproducible parameters"))
})

test_that("the planning prompt carries the injected skill block", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    # system_prompt is a formal of sclet_ai_call, not part of the ... the
    # sclet.ai.call mock receives, so capture at the call boundary itself.
    captured <- new.env(parent = emptyenv())
    fake_call <- function(task, context, ..., system_prompt = NULL) {
        captured$system_prompt <- system_prompt
        response <- list(
            answer = "mock plan",
            proposed_actions = list(list(id = "step1", action = "inspect_status")),
            success_criteria = list(),
            stop_conditions = list()
        )
        sclet:::sclet_ai_normalize_response(response,
            task = task, context = context,
            metadata = list(mock = TRUE, payload_policy = NULL))
    }
    testthat::with_mocked_bindings(
        {
            plan <- AIPlanAnalysis(sce, goal = "Cluster the cells", model = "mock")
            expect_true(!is.null(captured$system_prompt))
            expect_true(grepl("Single-cell analysis guidelines",
                captured$system_prompt, fixed = TRUE))
            expect_true(grepl("stage order", captured$system_prompt,
                fixed = TRUE))
            # the escape hatch really removes it, for cost-sensitive callers
            captured$system_prompt <- NULL
            old_opts <- options(sclet.ai.use_skills = FALSE)
            on.exit(options(old_opts), add = TRUE)
            invisible(AIPlanAnalysis(sce, goal = "Cluster the cells",
                model = "mock"))
            expect_false(grepl("Single-cell analysis guidelines",
                captured$system_prompt %||% "", fixed = TRUE))
        },
        sclet_ai_call = fake_call,
        .package = "sclet"
    )
})

# ---- M12-d: the four plan-quality dimensions and data_support must be
# ---- injected, not merely written into a file nobody reads.
#
# `options(sclet.ai.use_skills)` is a switch, and a skill edit that never
# reaches the prompt is invisible in every other test: the files parse, the
# selection matches, the prompt is built -- and the section the model needed
# is on disk only. These assertions pin the CONTENT reaching the model, which
# is what M12-d claims to deliver.

test_that("sclet-analysis carries the four plan-quality dimensions", {
    body <- ReadAISkill("sclet-analysis")$body
    expect_true(grepl("four testable dimensions", body, fixed = TRUE))
    for (dimension in c(
        "**Goal coverage.**",
        "**Data fitness.**",
        "**Round progression.**",
        "**Evidence increment.**"
    )) {
        expect_true(grepl(dimension, body, fixed = TRUE),
            info = paste("missing dimension:", dimension))
    }
    # the data-shape table is the actionable half of data fitness
    expect_true(grepl("Data shape decides the exploration order", body,
        fixed = TRUE))
    expect_true(grepl("`counts` plus `spliced`/`unspliced`", body,
        fixed = TRUE))
})

test_that("sclet-provenance carries how to read data_support", {
    body <- ReadAISkill("sclet-provenance")$body
    expect_true(grepl("analysis_story$data_support", body, fixed = TRUE))
    for (field in c("`shapes`", "`supports`", "`blocked`",
        "`unconfirmed_semantics`", "`unexplored`")) {
        expect_true(grepl(field, body, fixed = TRUE),
            info = paste("missing data_support field:", field))
    }
    # the one sentence that makes `unexplored` usable rather than a to-do list
    expect_true(grepl("data value not yet collected", body, fixed = TRUE))
})

test_that("sclet-driver states the repairs for the measured failure shapes", {
    body <- ReadAISkill("sclet-driver")$body
    # shape 1: an empty plan (1 of 9 Kang trials)
    expect_true(grepl("plan_empty", body, fixed = TRUE))
    expect_true(grepl("never the way to record", body, fixed = TRUE))
    # shape 2: a placeholder action used to drop a step (1 of 9)
    expect_true(grepl("action is not registered", body, fixed = TRUE))
    expect_true(grepl("\"not_run\"", body, fixed = TRUE))
    # shape 3: a dependency left pointing at a removed step (2 of 9)
    expect_true(grepl("depends on unknown step(s): step_N", body,
        fixed = TRUE))
    expect_true(grepl("exactly two honest repairs", body, fixed = TRUE))
})

test_that("data_support guidance is selected for planning goals, not only reporting ones", {
    # the regression this prevents: the provenance skill was selected by
    # report/comparison words, so a plain planning goal never received the
    # section that tells the model what the data can support next
    topics <- sclet_ai_goal_skill_topics("what else can this data support")
    expect_true("provenance" %in% topics)
    topics <- sclet_ai_goal_skill_topics("plan the next round")
    expect_true("provenance" %in% topics)
    # and the always-on analysis skill keeps the four dimensions available
    expect_true("analysis guidelines" %in%
        sclet_ai_goal_skill_topics("cluster the cells"))
})

test_that("sclet-driver states that a round ends with a rendered report", {
    body <- ReadAISkill("sclet-driver")$body
    expect_true(grepl("End every round by rendering the report", body, fixed = TRUE))
    expect_true(grepl("A round that is not", body, fixed = TRUE))
    expect_true(grepl("RenderAIReport", body, fixed = TRUE))
    # the next round's inputs must be named, not left to be guessed
    expect_true(grepl("data_support$unexplored", body, fixed = TRUE))
    expect_true(grepl("continue_blocked", body, fixed = TRUE))
})

test_that("sclet-analysis states where contrast levels come from", {
    body <- ReadAISkill("sclet-analysis")$body
    expect_true(grepl("contrast_levels", body, fixed = TRUE))
    expect_true(grepl("with `group = <column>`", body, fixed = TRUE))
    expect_true(grepl("identity_levels$codes", body, fixed = TRUE))
    # the recognisable symptom of inventing a level
    expect_true(grepl("you invented it", body, fixed = TRUE))
})

# ---- B20: the unexplored-response requirement must be in the prompt ---------
#
# Root cause of B20: `analysis_story$data_support$unexplored` was delivered in
# the context while not one line of the planning prompt mentioned that the block
# exists. Data without a covering instruction. Measured: 4 of 6 trials on the
# ambient fixture ignored the suggestion and never mentioned it. The fix is a
# RESPONSE requirement -- each entry planned or declined by name -- and these
# pins hold that wording in both prompts, because a prompt edit that never
# reaches the model is invisible to every other test.

test_that("the planning prompt requires each unexplored entry to be planned or declined", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 5L, ncol = 4L))
    )
    captured <- new.env(parent = emptyenv())
    fake_call <- function(task, context, ..., system_prompt = NULL) {
        captured$system_prompt <- system_prompt
        response <- list(
            answer = "mock",
            proposed_actions = list(list(id = "step1", action = "inspect_status")),
            success_criteria = list(),
            stop_conditions = list()
        )
        sclet:::sclet_ai_normalize_response(response, task = task,
            context = context, metadata = list(mock = TRUE))
    }
    testthat::with_mocked_bindings(
        {
            invisible(AIPlanAnalysis(sce, goal = "cluster the cells", model = "mock"))
            prompt <- captured$system_prompt %||% ""
            expect_true(grepl("unexplored", prompt, fixed = TRUE))
            # the response requirement, not a "you may find this useful"
            expect_true(grepl("EACH entry must be either PLANNED", prompt,
                fixed = TRUE))
            expect_true(grepl("explicitly DECLINED", prompt, fixed = TRUE))
            # the other two blocks are named so the model knows where to look
            expect_true(grepl("`blocked`", prompt, fixed = TRUE))
            expect_true(grepl("`shapes`", prompt, fixed = TRUE))
        },
        sclet_ai_call = fake_call,
        .package = "sclet"
    )
})

test_that("the tool-loop prompt carries the same requirement", {
    line <- sclet:::sclet_ai_unexplored_prompt_line()
    expect_true(grepl("EACH entry must be either PLANNED", line, fixed = TRUE))
    expect_true(grepl("treated as an oversight", line, fixed = TRUE))
    # declining is legitimate: the requirement must not read as "always plan it"
    expect_true(grepl("Declining is legitimate", line, fixed = TRUE))
})
