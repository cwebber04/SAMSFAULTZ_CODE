#resources
http://gogojimmy.net/2012/01/17/how-to-use-git-1-git-basic/
#reveal hidden files:
defaults write com.apple.finder AppleShowAllFiles TRUE;\killall Finder
#hide hidden files
defaults write com.apple.finder AppleShowAllFiles FALSE;\killall Finder
#configurization
git config --list
git config --global user.name "TimothyTLee"
git config --global user.email "webber04@gmail.com"
#alias function:
git config --global alias.st status
git config --global alias.co checkout
git config --global alias.cm commit
git config --global alias.br branch
#ignore space:
git config --global apply.whitespace nowarn
#color interface:
git config --global color.ui true
#ignore tracking
vi .gitignore
# common command
git log
git fetch
git pull
git push
#repository
git init
git remote add origin https://github.com/cwebber04/SAMSFAULTZ_CODE.git
git push -u origin master
git clone https://github.com/cwebber04/Obspy.git Obspy #duplicate on github
git log [--stat]/[-a]
#define the current verison of the code
git branch git
git checkout git
#check if new items in the depositonary
git status
git add "..."
git add . #add all untracked files
git all -i #add by interactive mode
git status
git commit
git commit -m "..." #quick commit message
#organize all branches
git rebase master
git diff master git
#merge and push the newer version to github
git branch
git checkout master
git merge git
git push
#cancel merge
git reset --hard ORIG_HEAD
#cancel add
git reset HEAD git_basics.csh
#recover file from last commit
git checkout -- git_basics.csh
#revise commit message
git commit --amend
