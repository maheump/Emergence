# Finger tracking analyses

The script ```Emergence_FTA_PreprocData``` preprocesses the raw data, it has to be run only once. It generates a MATLAB file with data from all the subjects and corresponding inference from the full ideal observer of the task.

The script ```Emergence_FTA_LoadData``` loads the preprocessed data in the MATLAB workspace.

All other scripts are analyses scripts that can be run either after calling ```Emergence_FTA_LoadData``` as follows.

```
>> Emergence_FTA_LoadData
Loading data... Done! Bad subjects are excluded.
>> Emergence_FTA_HypothesisWeighting1
```

The script ```Emergence_FTA_RunPipeline``` runs the entire analysis pipeline.
