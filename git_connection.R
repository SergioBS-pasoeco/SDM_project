
library(gert)
library(usethis)

create_github_token()


# 1) Ensure identity (skip if already set)
usethis::use_git_config(user.name = "SergioBS-pasoeco", user.email = "sergiobolivars@gmail.com")
