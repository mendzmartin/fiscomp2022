# Outputs

+ prob01_01.o
```
...$ gfortran -g3 -Wall -Wextra -Wconversion -o prob01_01 prob01.f90 && ./prob01_01
prob01.f90:12:16:

   12 |     res = 37/3      ! prob01_01.o
      |                1
Warning: Nonconforming tab character at (1) [-Wtabs]
prob01.f90:12:12:

   12 |     res = 37/3      ! prob01_01.o
      |            1
Warning: Integer division truncated to constant ‘12’ at (1) [-Winteger-division]
 result =    12.000000000000000
```
+ prob01_02.o
```
...$ gfortran -g3 -Wall -Wextra -Wconversion prob01.f90 -o prob01_02 && ./prob01_02
prob01.f90:13:2:

   13 |  res = 37 + 17 / 3    ! prob01_02.o
      |  1
Warning: Nonconforming tab character at (1) [-Wtabs]
prob01.f90:13:14:

   13 |  res = 37 + 17 / 3    ! prob01_02.o
      |              1
Warning: Integer division truncated to constant ‘5’ at (1) [-Winteger-division]
 result =    42.000000000000000
```
+ prob01_03.o
```
gfortran -g3 -Wall -Wextra -Wconversion prob01.f90 -o prob01_03 && ./prob01_03
prob01.f90:14:2:

   14 |  res = 28 / 3 / 4    ! prob01_03.o
      |  1
Warning: Nonconforming tab character at (1) [-Wtabs]
prob01.f90:14:9:

   14 |  res = 28 / 3 / 4    ! prob01_03.o
      |         1
Warning: Integer division truncated to constant ‘9’ at (1) [-Winteger-division]
prob01.f90:14:11:

   14 |  res = 28 / 3 / 4    ! prob01_03.o
      |           1
Warning: Integer division truncated to constant ‘2’ at (1) [-Winteger-division]
 result =    2.0000000000000000
```
+ prob01_04.o
