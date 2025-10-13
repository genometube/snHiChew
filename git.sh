######## ssh connection
ssh-keygen -t ed25519 -C "genometube2@gmail.com"
# eval "$(ssh-agent -s)"
# ssh-add ~/.ssh/id_ed25519

#Copy Public Key to Clipboard:
# cat ~/.ssh/id_ed25519.pub   
cat /root/.ssh/id_ed25519.pub

# Add SSH Key to GitHub:

# Go to GitHub Settings > SSH and GPG keys
# Click "New SSH key"
# Paste your public key and save

REPO_URL="git@github.com:genometube/MICC-seq.git"
git clone $REPO_URL

# Test SSH Connection
ssh -T git@github.com

git remote set-url origin git@github.com:genometube/snHiChew.git

# Navigate to your repo
cd snHiChew

# git=/mnt/software/anaconda3/envs/R4_4/bin/git
git=/usr/bin/git
$git remote -v 
$git add --all .
$git commit -m "update"
$git push origin main

# sync remote to local
git pull origin main

# # new section
# mkdir test_section
# touch test_section/test.md
# git add test_section
# git commit -m "Create new subfolder and file"
# git push origin master

# Your public key has been saved in /root/.ssh/id_ed25519.pub.
# The key fingerprint is:
# SHA256:RK1hPAsE2TBhZgYd6nggdHEngC6PXf3J6J2EDhS/nbk genometube2@gmail.com
# The key's randomart image is:
# +--[ED25519 256]--+
# | .o*&@.o..       |
# |...*=.=.= .      |
# |+ .  + o.=       |
# |o=  o o.o        |
# |o+oo   BS+       |
# |..o . + O        |
# |     + o o       |
# |      o E        |
# |                 |
# +----[SHA256]-----+