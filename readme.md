# Supporting Data
This repo contains supporting data for papers published from the Holehouse lab. In addition to sharing our data here, below we provide some brief instructions on setting up your own supporting data GitHub repository. 


## Setting up your own  supporting data repository
###### Last updated May 2021
The following instructions provide a quick summary on how to set up your own supporting data repository. In general GitHub repos suffer from the downside of being non DOIed and having a max size of ~1GB. If having a DOI for your data is essential, and/or you need more space than 1 GB, then [Zenodo](https://zenodo.org/) is an alternative, complementary, and excellent option (and by many people might be recommended). 

### Why use `git` over Zenodo?
Using git makes it easy to quickly update/add/refine data, which may be a benefit (can be responsive to questions from the community, add additional info in a public way, clarify things) BUT may also be a downside (changes from the original published supporting data can happen). 

If the goal is the 'snapshot' the info when published Zenodo might be better. 

If the goal is to provide a set of data/analysis which may evolve and change over time, git may be better.

If the goal is to share code then GitHub is almost certainly better.

The BEST case scenario is probably to submit data for a paper via a Zenodo repo and also have an accessible and evolving dataset on GitHub... But, that might be overkill. 

### Setting up a GitHub supporting data repository

The instructions below are valid of *nix systems (macOS or Linux/unix). They assume familiarity with the command-line/terminal and are pretty terse. The instructions also try to avoid getting into `git` as a technology as best is possible...  

1. If you don't yet have one, create a GitHub account!

2. Create a new repository in your account - [see instructions here](https://docs.github.com/en/github/getting-started-with-github/create-a-repo). Do not include any additional info (e.g. license files etc) here.

3. Check if `git` is installed on your local machine (it really should be). To do this, open the terminal and run:

		git --version
		
	This should return information on the current version of git installed. If you git a `command not found` error you'll need to install git.
	
4. Next you need to link your local computer to your GitHub account. See [the information here](https://docs.github.com/en/github/authenticating-to-github/about-authentication-to-github) on authenticating, which can be done by SSH-keys or various.

5. Finally, choose a location on your local computer and from the terminal run
 
		git clone git@github.com:<username>/<repo name>.git
	This info is copy-able from the big green CODE button in the top right of the repository you just created, e.g. for our repo this is.
	
		git clone git@github.com:holehouse-lab/supportingdata.git


	Once you've cloned in the repo this creates a directory on your local computer that is linked to the GitHub repository. 

6. 				

	The remaining workflow involves copying files into this directory, **staging** them in the local repository using `git add`, **committing** staged files via `git commit` and then **pushing** the local repository to github using `git push`.
	
	##### Adding files to the staging area
	`git add` lets you tell your repository that there are files you want to upload. Specifically, once a file is copied into the git repository directory, it can be added using:
				
		git add <filename>
		
		# You can also select all files recursively in the current dir via
		git add .
		
	As you add files they go into a staging area - basically a place where they hang out before you `commit` (which is the act of updating the git history). You can unstage files if you want via:
	
		git restore --staged <filename>

	Importantly, just copying a file into the directory DOES NOT automatically add it. You can see which (new) files have or have no been added/committed by running:
	
		git status
		
	When you add a file, it takes a snapshot of the file _at the moment it was *added* so if you make further changes you have to re-add. The staging area means in a complicated project with many files you can be selective about which ones will get sent up to the GitHub repo.
	
	Once you're happy with the files you can commit then:
	
		git commit -m "Commit message"
		
	The commit message is a short description of what was changed/included. Committing is the act of saying "Yes, these files should be added to the version history of the file". 
	
	Finally, once committed you can run:
	
		git push
		
	This transmits the info from your local machine to the remote repository (i.e. on GitHub) and your files should now be visible online. 
	
	This workflow lets you update/add new files at your convenience, and all from the comfort of the commandline.
	
	
	
### Epilog
**git** is an amazing technology. I highly recommend [this video introduction](https://missing.csail.mit.edu/2020/version-control/) to git if you're unfamiliar, as it helps explain the conceptual basis which makes really understand what we're doing much, much easier.
			
			