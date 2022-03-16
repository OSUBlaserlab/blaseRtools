
# run these commands to rewind to a prior "good" commit ----------------------

# make sure git status is "clean" (all changes committed) before rewinding
# gert::git_log() #find the id of the good commit
# blaseRtemplates::git_rewind_to(commit = "<good commit id>")



# # Run these commands regularly for branching, updating and merging --------

# For more complicated situations, you should
# consider using the terminal or a gui program.

# create a working branch for your day's work
blaseRtemplates::git_easy_branch(branch = "brad_working")

# save, add and commit your work but don't push
gert::git_add("*")
gert::git_commit("")

# frequently update your working branch from main or master branch
# this will first update main or master from remote
blaseRtemplates::git_update_branch()


# once you are done with your day's work, merge back into main
blaseRtemplates::git_safe_merge()

# remember to delete your branch when you are done merging:
gert::git_branch_delete(branch = "brad_working")

# remember to push your changes to github so we can all get them:
blaseRtemplates::git_push_all()



# conflict resolution -----------------------------------------------------

# # any conflicting updates will be marked and the files will need to be edited
# # to resolve the conflicts.
#
# # in extreme circumstances you may need to accept all changes in one file from one source
# # to accept all changes from the main or master branch (pulled from remote) run in the terminal:
# # git checkout --ours
#
# # to accept all changes from working branch
# # git checkout --theirs
# #
# # then run in the terminal
# # git add .
# # git rebase --continue
#
# # you will be presented with a commit message. Enter :wq in the terminal to save if using Vim.





