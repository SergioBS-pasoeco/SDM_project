library(gert)

gert::git_remote_list() # the origin is the SDM_project.git
gert::git_fetch()
gert::git_push(set_upstream = TRUE)


repo_path <- git_find()  # errors if not in a git repo

git_info()               # current branch, head, repo path
git_status()             # staged/unstaged changes
git_branch_list()        # local branches
git_remote_list()        # remotes (origin, etc.)
git_log(max = 10)        # recent commits
git_config()

git_find()

gert::git_add("scripts/read_occurrences.R")
gert::git_commit("Add script to read and standardize occurrence data")
usethis::use_github()
# First push from this branch:
gert::git_push(set_upstream = TRUE)



# Create and switch to a new branch
gert::git_branch_create("add-read-occurrences")

# Push the new branch to GitHub
gert::git_push(set_upstream = TRUE)
