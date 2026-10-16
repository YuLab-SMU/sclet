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
