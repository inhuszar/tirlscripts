# tirlscripts-oxford-mnd

Scripts in this repository constitute the 4-stage histology-to-MRI registration 
pipeline that was originally introduced by _Huszar et al_(NeuroImage, 2023) for 
a post-mortem imaging dataset of a cohort of Motor Neuron Disease (MND) 
patients and disease-free controls.

-------------

If you use any of these tools, please cite the following peer-reviewed journal 
article:

**IN Huszar**, M Pallebage-Gamarallage, S Bangerter-Christensen, H Brooks, 
S Fitzgibbon, S Foxley, M Hiemstra, AFD Howard, S Jbabdi, DZL Kor, A Leonte, 
J Mollink, A Smart, BC Tendler, MR Turner, O Ansorge, KL Miller, 
and M Jenkinson: _Tensor image registration library: Deformable registration of 
stand‐alone histology images to whole‐brain post‐mortem MRI data._
NeuroImage, Volume 265, 2023, 119792, ISSN 1053-8119, DOI: 
<a href=https://dx.doi.org/10.1016/j.neuroimage.2022.119792>10.1016/j.neuroimage.2022.119792</a>.

--------------

When running the pipeline, the expectation is that you will not need to change 
anything in the registration scripts. However, you will need to edit the YAML 
configuration file (*.yml) that accompanies each of the registration scripts, 
to specify the input files and the output folder. If the inputs are 
considerably different from the images used in our study, you may also need to 
adapt some of the numerical registration parameters. The default values found 
in the configuration files represent empirical optima for our dataset. Further 
information on the registration parameters can be found within the comments of 
the configuration files.

A registration script with a configuration file should be executed from the 
Terminal as follows (make sure that the tirl environment is active):

```bash
tirl mnd.h2h --config histology_to_histology.yml --verbose
```

Note that `mnd.h2h` is an alias for the histology_to_histology.py script, which 
is installed by this module. Installed scripts and their aliases are listed in 
`~/.config/tirlscripts.yml`. By editing this file, you can change the aliases 
of the installed scripts, or add locally modified versions of the scripts under 
a custom alias. Configuration files are not aliased, therefore these must be 
referenced by their full path in the command.

## The Pipeline

### Stage 0: histology_to_histology.py

- Registers different histology stains of the same tissue section 
(or a close proxy).


### Stage 1: histology_to_block.py

- Registers a histology section to a blockface photograph.


### Stage 2: block_to_slice.py

- Registers a blockface photograph to a brain slice photograph. The sampling 
region must be specified manually or using the output of the next script.

- **find_sites.py** : Based on a sequence of brain slice photographs, 
identifies regions where tissue blocks were removed from the brain slice, i.e. insertion 
sites for the blockface photographs. 


### Stage 3: slice_to_volume.py

- Registers an individual brain slice photograph to a 3D MRI volume in a direct, 
2D-to-3D scheme.


### Stage 4: histology_to_volume.py

- Based on the optimised transformations from Stage 1-3, this script fine-tunes 
the alignment of the histology image in 3D MRI space. 


## Need help?

For further information and guidance on running the pipeline, please consult 
the tutorials in the `tirl-tutorials` repository, or send me an email: 
istvan.huszar at ndcn.ox.ac.uk.
