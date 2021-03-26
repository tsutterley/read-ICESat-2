convert_delta_time.py
=====================

 - Converts time from delta seconds into Julian and year-decimal

#### Calling Sequence
```python
from icesat2_toolkit.convert_delta_time import convert_delta_time
t_date = convert_delta_time(delta_time)['decimal']
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/convert_delta_time.py)

#### Arguments
 1. `delta_time`: seconds since gps_epoch

#### Keyword arguments
 - `gps_epoch`: seconds between delta_time and GPS epoch (1980-01-06T00:00:00)

#### Returns
 - `julian`: time in Julian days
 - `decimal`: time in year-decimal
