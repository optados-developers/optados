module od_build
 implicit none
 
 private ! unless otherwise stated
 
 type, public ::  build_info_type
  character(len=20) :: build=''
  character(len=20) :: compiler='gfortran'
  character(len=20) :: build_type='debug'
  character(len=20) :: comms_arch='serial'
  character(len=20) :: source_time='17:36:13'
  character(len=20) :: source_date='Mon 19 Dec 2022'
  character(len=20) :: compile_date='Mon 19 Dec 2022'
  character(len=20) :: compile_time='17:36 GMT'
 end type build_info_type
 type(build_info_type), public, save :: build_info
endmodule od_build
