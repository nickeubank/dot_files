

# Git completion
if [ -f ~/.git-completion.bash ]; then
   source ~/.git-completion.bash
   source ~/.git-prompt.sh
   export PS1='\W$(__git_ps1 "(%s)")\$'
fi


# added by Anaconda3 2.3.0 installer
export PATH="/Users/Nick/anaconda/bin:$PATH"

alias pdev="cd ~/github/pandas; source activate pandas_dev"
alias pdocs="cd ~/github/pandas/doc; echo python make.py --single NAME"

alias pupdate="git fetch upstream;git checkout master;git merge upstream/master;git branch -a"
alias polgeo="cd ~/dropbox/gis_in_r;"

alias pss="cd ~/github/python-for-social-scientists"

##
# Your previous /Users/Nick/.bash_profile file was backed up as /Users/Nick/.bash_profile.macports-saved_2015-11-04_at_13:42:48
##

# MacPorts Installer addition on 2015-11-04_at_13:42:48: adding an appropriate PATH variable for use with MacPorts.
export PATH="/opt/local/bin:/opt/local/sbin:$PATH"
# Finished adapting your PATH environment variable for use with MacPorts.

