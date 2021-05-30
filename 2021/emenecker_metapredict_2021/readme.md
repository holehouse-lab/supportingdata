## Additional Supporting Information
###### Last updated 2021-05-29

This directory contains information for the manuscript

**metapredict: a fast, accurate, and easy-to-use cross-platform predictor of consensus disorder** by Ryan J. Emenecker, Daniel Griffith, and Alex S. Holehouse.

### Useful links:

**Full documentation for metapredict:** [https://metapredict.readthedocs.io/](https://metapredict.readthedocs.io/)

**The metapredict GitHub repository:** [https://github.com/idptools/metapredict](https://github.com/idptools/metapredict)

**The metapredict webserver:** [http://metapredict.net/](http://metapredict.net/).

## Contents

### `/performance`
Contains code for assessing `metapredicts` computational performance (i.e. residues-per-second). Note this requires the package [protfasta](https://protfasta.readthedocs.io/) to run (as does metapredict).

### `/demo`
This contains a demo-notebook with a subset of the analysis that the metapredict Python package offers. The contents of this directory can be downloaded and run locally.

### `/comparison_figure_3`
Code to generate the Figure 3 in the metapredict paper

### `/data`
All of the CAID and CheZOD data, both scores and sequences, for easy future comparison. If you use either of these you MUST cite [1] and [2]. We are providing these files to simpify future analysis, but they originate from these papers and should be cited as such.


## FAQs
*Written by Alex*

#### I can't get Metapredict installed  - what do I do!?
Firstly, make sure you can install other Python packages. We _highly_ recommend installing **metapredict** into a `conda` environment. If that idea is confusing to you, please check out how to get [conda up and running](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) and perhaps check out [a tutorial](https://www.youtube.com/watch?v=1VVCd0eSkYc) or [two](https://www.youtube.com/watch?v=qn5zfdJtcYc).

Assuming you have `conda` set up and an environment loaded:

	pip install metapredict
	
should _just work_. 

If you can install other packages but not **metapredict**, please either contact [Alex or Ryan directly](https://www.holehouselab.com/team) OR [raise a GitHub Issue](https://github.com/idptools/metapredict/issues) on the **metapredict** GitHub page! If you're struggling with getting `conda` up and running in general, please get that sorted first _before_ you contact (although, if you're really stuck, Alex can probably help you because he's a sucker).

#### I found a bug (I think) in metapredict!
Please either contact [Alex or Ryan directly](https://www.holehouselab.com/team) OR [raise a GitHub Issue](https://github.com/idptools/metapredict/issues) on the **metapredict** GitHub page! Don't hesitate if you found something that seems bug-like; worst case our docs are confusing/wrong and we should fix that. "Best" case, you found a bug and we need to fix it! Importantly, if you found something confusing _we_ need to make it clearer - *i.e.,* the burden of user experience is on us (the developers), not you (the users)!

#### Why would you make ANOTHER disorder predictor? Seriously.
The honest answer is we needed a fast, open source, easy to deploy/use and cross-platform disorder predictor that worked in Python, and we couldn't find one that met these requirements. This was the motivation for **metapredict**. If OTHER people find it useful, that's great, but, at some basic level we built **metapredict** for ourselves and found it so useful we wanted to share it with the world. If you don't like it, don't use it!

#### You did not cite my predictor and it's better than metapredict!
We're sorry! As mentioned, we built **metapredict** to enable us to have a fast, accurate, and portable predictor of disorder to get *other* stuff done. We did not build it to poke people in the eye or imply other predictors are bad or that we are the best. We *do* think **metapredict** makes it _easy_ to integrate disorder prediction into whatever workflow you're doing, and, in our experience, this was a set of tools missing from the bioinformatics ecosystem. However, again, please don't be offended if we failed to cite/use/compare your predictor.

#### I want a new feature in metapredict!
Please - feel free to add a feature request via [GitHub Issue](https://github.com/idptools/metapredict/issues). We'd actually love this - we really want **metapredict** to make people's lives easier, and we built **metapredict** to make it easy to add features. This is true for all three "branches" of **metapredict** (Python package, commandline tools, web server).

You can also branch, add a new feature, and issue a pull request. If you do this, please make sure (1) you include tests for your new feature (testing is done via `pytest` and (2) you include documentation, and the updated docs build correctly in Sphinx.

## References
[1] Necci, M., D. Piovesan, CAID Predictors, DisProt Curators, and S.C.E. Tosatto. 2021. Critical assessment of protein intrinsic disorder prediction. Nat. Methods.

[2] Nielsen, J.T., and F.A.A. Mulder. 2019. Quality and bias of protein disorder predictors. Sci. Rep. 9:1–11.


