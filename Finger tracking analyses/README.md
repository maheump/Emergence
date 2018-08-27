# Finger tracking analyses

The script ```Emergence_FTA_PreprocData``` preprocesses the raw data, it has to be run only once. It generates a MATLAB file with data from all the subjects and corresponding inference from the full ideal observer of the task.

The script ```Emergence_FTA_LoadData``` loads the preprocessed data in the MATLAB workspace.

The script ```Emergence_FTA_RunPipeline``` runs the entire analysis pipeline.

All other scripts are analyses scripts that can be run either after calling ```Emergence_FTA_LoadData``` or by running the entire pipeline script ```Emergence_FTA_RunPipeline```.
