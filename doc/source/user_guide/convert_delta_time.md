convert_delta_time.py
=====================

 - Converts time from delta seconds into Julian and year-decimal  

#### Calling Sequence
```python
from icesat2_toolkit.convert_delta_time import convert_delta_time
t_date = convert_delta_time(delta_time)['decimal']
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/master/icesat2_toolkit/convert_delta_time.py)

#### Inputs
 1. `delta_time`: seconds since gps_epoch  

#### Options
 - `gps_epoch`: seconds between delta_time and GPS epoch (1980-01-06T00:00:00)  

#### Outputs
 - `julian`: time in Julian days  
 - `decimal`: time in year-decimal  
