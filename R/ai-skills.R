# AI-readable analysis skills, shipped with the package.
#
# WHY THESE LIVE IN inst/skills
#
# The 2026-10-16 real-provider benchmark sweep (sclet-benchmark-ai) found
# that the model's failures were not engine failures -- the engine validated,
# refused and recorded correctly -- but knowledge failures: it ordered
# run_rare_cell_detection before its evidence producers 2/2 times, it passed
# an OrgDb where an ontology database belongs, it reached for the reference
# syntax it had been taught only after being corrected. That knowledge exists
# in the single-cell literature as guidelines (Heumos 2023 et al.) and in
# provenance and claim-calibration standards, but none of it was in the
# planning context, so the model had no path to knowing it.
#
# These files are that knowledge, in the SKILL.md format agents already read
# (YAML frontmatter with name + description, then body), so they serve two
# consumers from one source:
#
# - an external AI agent (Claude Code, or any skill-directory reader) can
#   load them directly from the installed package's
#   system.file("skills", package = "sclet") directory;
# - sclet's own planning and interpretation prompts can inject the body of
#   the skills relevant to the task at hand, through
#   sclet_ai_skill_context().
#
# One source of truth; two injection paths. The alternative -- prompting
# rules inline in R strings -- is how the current prompt got its parameter
# rules, and it cannot carry the volume of domain knowledge these files do
# without bloating every call.

#' List the analysis skills shipped with sclet
#'
#' Skills are plain SKILL.md files under `inst/skills/<name>/`, readable by
#' any AI agent, and by the internal \code{sclet_ai_skill_context()}.
#'
#' @return A character vector of skill directory names, invisibly.
#' @examples
#' skills <- AISkills()
#' length(skills) > 0
#' @export
AISkills <- function() {
    root <- system.file("skills", package = "sclet")
    dirs <- list.dirs(root, recursive = FALSE)
    sort(basename(dirs))
}

#' Read one shipped skill as a structured list
#'
#' Parses the YAML frontmatter (name, description) and returns the body
#' separately, so a caller can show a model the description cheaply and the
#' body only when the skill is relevant.
#'
#' @param skill Skill name, as returned by [AISkills()].
#' @return List with `name`, `description`, `body` and `path`.
#' @examples
#' skill <- ReadAISkill(AISkills()[[1]])
#' skill$name
#' nchar(skill$body) > 0
#' @export
ReadAISkill <- function(skill) {
    stopifnot(is.character(skill), length(skill) == 1L, nzchar(skill))
    path <- system.file("skills", skill, "SKILL.md", package = "sclet")
    if (!nzchar(path)) {
        stop("No such skill: '", skill, "'. Available: ",
            paste(AISkills(), collapse = ", "), call. = FALSE)
    }
    lines <- readLines(path, warn = FALSE)
    fm_end <- if (length(lines) && identical(trimws(lines[[1L]]), "---")) {
        end <- which(trimws(lines) == "---")[[2L]]
        end
    } else {
        0L
    }
    frontmatter <- if (fm_end >= 2L) lines[2L:(fm_end - 1L)] else character()
    body <- if (length(lines) > fm_end) lines[(fm_end + 1L):length(lines)] else character()
    get_field <- function(key) {
        start <- which(grepl(paste0("^", key, ":"), frontmatter))
        if (!length(start)) return("")
        start <- start[[1L]]
        first <- sub(paste0("^", key, ":[[:space:]]*"), "", frontmatter[[start]])
        first <- sub("[[:space:]]+$", "", first)
        # A folded block scalar (`description: >`) carries the text on the
        # indented continuation lines, so a first line holding just the
        # marker means "keep reading". A value on the marker line itself is
        # the whole value.
        if (nzchar(first) && !first %in% c(">", "|", ">-", "|-")) {
            return(trimws(first))
        }
        rest <- if (start < length(frontmatter)) {
            frontmatter[seq.int(start + 1L, length(frontmatter))]
        } else {
            character()
        }
        keep <- rest[grepl("^[[:space:]]+", rest) | nzchar(trimws(rest)) &
            grepl("^[[:space:]]", rest)]
        # Stop at the next field key, not just blank lines: frontmatter is
        # short, but a description that mentions "name:" must not swallow it.
        keep <- keep[!grepl("^[[:space:]]*[a-z_]+:", keep)]
        trimws(paste(trimws(keep), collapse = " "))
    }
    list(
        name = get_field("name") %||% skill,
        description = get_field("description"),
        body = paste(body, collapse = "\n"),
        path = path
    )
}

#' Build the skill context block for an AI prompt
#'
#' Selection is by substring match between `topics` and each skill's
#' frontmatter description, so a planning call for an enrichment-chain goal
#' can pull in the analysis skill without the claims skill. When every skill
#' matches (or `topics` is NULL), all bodies are included -- cheap enough
#' that precision of selection matters less than the injection itself.
#'
#' @param topics Character vector of topic keywords, matched against skill
#'   descriptions (case-insensitive). NULL includes every skill.
#' @param format `"bodies"` returns full skill bodies ready to inject into a
#'   system prompt; `"index"` returns one description line per skill, for a
#'   context where the model may ask for a body by name.
#' @return A single string, "" when nothing matches.
#' @keywords internal
sclet_ai_skill_context <- function(topics = NULL, format = c("bodies", "index")) {
    format <- match.arg(format)
    skills <- AISkills()
    if (!length(skills)) return("")
    parsed <- lapply(skills, ReadAISkill)
    names(parsed) <- skills
    if (!is.null(topics)) {
        topics <- tolower(as.character(topics))
        hit <- vapply(parsed, function(s) {
            text <- tolower(paste(s$name, s$description))
            any(vapply(topics, function(t) grepl(t, text, fixed = TRUE), logical(1)))
        }, logical(1))
        parsed <- parsed[hit]
    }
    if (!length(parsed)) return("")
    if (identical(format, "index")) {
        paste(vapply(parsed, function(s) {
            paste0("- ", s$name, ": ", s$description)
        }, character(1)), collapse = "\n")
    } else {
        paste(vapply(parsed, function(s) {
            paste0("## Skill: ", s$name, "\n\n", s$body)
        }, character(1)), collapse = "\n\n")
    }
}

#' Map an analysis goal to the skills worth injecting
#'
#' Keyword-based and deliberately over-inclusive: including a skill that was
#' not needed costs characters, while excluding the one that was needed costs
#' a repair round. The topics are matched against skill descriptions, so a
#' new skill participates in selection as soon as its description mentions
#' the right words.
#'
#' @param goal The user's analysis goal string.
#' @return Character vector of topics for [sclet_ai_skill_context()], or NULL
#'   to include everything.
sclet_ai_goal_skill_topics <- function(goal) {
    goal <- tolower(paste(goal, collapse = " "))
    topics <- character()
    # analysis skill: always relevant to planning -- it owns stage order,
    # action choice and parameter discipline, which every plan touches.
    topics <- c(topics, "analysis guidelines")
    # provenance skill: relevant when the plan reports, compares or
    # reproduces, and for criteria-bearing plans.
    #
    # `plan`, `next round`, `unexplored` and `data support` were added for M12:
    # the provenance skill now also carries how to read
    # `analysis_story$data_support`, which is the input to choosing what the
    # NEXT round should do. A goal that never says "report" still needs it --
    # "what else can this data support" is a planning question, not a
    # reporting one.
    if (grepl("report|reproduc|provenance|compar|criterion|success|version|seed|plan|next round|unexplored|data support",
        goal)) {
        topics <- c(topics, "provenance")
    }
    # claims skill: relevant to interpretive calls and to any plan whose
    # output is biological language rather than artifacts.
    if (grepl("marker|annotation|cell type|interpre|claim|evidence|traject|fate|signall|communicat|enrich|finding|report|rare",
        goal)) {
        topics <- c(topics, "claim")
    }
    if (!length(topics)) NULL else topics
}
