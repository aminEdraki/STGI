# STGI
This is a MATLAB implementation of STGI, the speech intelligibility prediction algorithm proposed in [1]. The algorithm takes the time-aligned clean and degraded/processed speech signals as inputs. The output ```rho``` is expected to have a monotonically increasing relation with speech intelligibility. Please refer to Sec. V of [2] cited below for guidance on interpreting algorithm output.

## Installation
Before using the STGI function, the ```utils``` folder has to be added to the search path for the current MATLAB session.


## Usage
The function ```stgi``` takes three inputs:
```MATLAB
rho = stgi(clean_speech, degraded_speech, Fs);
```

* ```clean_speech```: An array containing a single-channel clean (reference) speech signal.
* ```degraded_speech```: An array containing a single-channel degraded/processed speech signal.
* ```Fs```: The sampling frequency of the input signals in ```Hz```.

Note that the clean and degraded speech signals must be time-aligned and of the same length.


## References
If you use STGI, please cite the reference [1] below:
```
[1] A. Edraki, W.-Y. Chan, J. Jensen, & D. Fogerty, “A Spectro-Temporal Glimpsing Index (STGI) for Speech Intelligibility Prediction," Proc. Interspeech, 5 pages, Aug 2021.
[2] A. Edraki, W.-Y. Chan, J. Jensen, & D. Fogerty, “Speech Intelligibility Prediction Using Spectro-Temporal Modulation Analysis,” IEEE/ACM Trans. Audio, Speech, & Language Processing, vol. 29, pp. 210-225, 2021.
```

## License
[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)
