* New branches cannot be made on MFlowCode/MFC, they are made on forks
* PRs:
    * made using AI tools like Claude Code and Codex should say so.
    * are made from those MFC forks
    * that change CFD results need verification that the PR is correct
    * follow template
    * that break a feature but promise a followup PR to fix it are rejected
* Commands:
    * MFC should almost always build and run using the ./mfc.sh command
    * Running mfc.sh commands can create a lock file in build/ that is sticky, be careful
* Programming and Design:
    * Keep this research fork aligned with MFlowCode/MFC master refactors so local features remain practical to submit upstream.
    * Extend upstream module boundaries, EOS and parameter registries, code generators, and GPU macros. Keep shared-core changes minimal and custom physics in focused modules.
    * Fetch the upstream remote before integration, preserve upstream merge ancestry and local physics, and validate each change without regenerating unexplained golden differences.
    * Keep feature changes separable into small, independently validated upstream PRs.
    * New code should follow the DRY principle and also make side-effect code DRY as well
    * Comments should be as short as possible without sacrificing value
    * GPU macros should follow the source's existing GPU macro principles and patterns
    * Functions/subroutines/modules shorter is better while being correct, fast, and separating concerns
