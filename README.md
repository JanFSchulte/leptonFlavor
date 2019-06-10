## Measurement of lepton flavor universality
- The plotLepFlavor.py is used to plot the LFU for DY and Other backgrounds (by changing background setting in main function), cmd: `python plotLeptFlavor.py -r --inverseAE`
- The plotJetsFlavor.py is used to plot Jets LFU, due to laziness in adding if-statements, cmd: `python plotJetsFlavor.py -r --inverseAE`.
- plotAllMCFlavor.py is used to compare MC to Data, cmd: `python plotAllMCFlavor.py -r -d --inverseAE`.

## Set limits
- setlimits.py is used to set a limit on LFU with signal=MC, bkg=Data.
