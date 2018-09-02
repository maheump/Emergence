# Finger tracking analyses

The script ```Emergence_FTA_PreprocData``` preprocesses the raw data, it has to be run only once. It generates a MATLAB file with data from all the subjects and corresponding inference from the full ideal observer of the task.

The script ```Emergence_FTA_LoadData``` loads the preprocessed data in the MATLAB workspace.

All other scripts are analyses scripts that can be run after calling ```Emergence_FTA_LoadData``` as follows.

```
>> Emergence_FTA_LoadData
Loading data... Done! Bad subjects are excluded.
>> Emergence_FTA_DistributionOfBeliefs
```

The script ```Emergence_FTA_RunPipeline``` runs the entire analysis pipeline (all the analyses scripts at once).

Note that some of the scripts can be run either on the actual behavior (default) or on the inference from the ideal observer by switching the ```D``` variable the following way.

```
>> D = G; % data from the subjects (default)
>> Emergence_FTA_DistributionOfBeliefs % distribution of finger positions (from the subjects)
>> D = IO; % data from the ideal observer
>> Emergence_FTA_DistributionOfBeliefs % distribution of posterior beliefs (from the IO)
```

Please refer to ```Emergence_FTA_RunPipeline``` in order to identify which scripts can be run based either on the subjects or the ideal observer.
