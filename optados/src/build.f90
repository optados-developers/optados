module od_build
 implicit none
 
 private ! unless otherwise stated
 
 type, public ::  build_info_type
  character(len=20) :: build=''
  character(len=20) :: compiler='gfortran'
  character(len=20) :: build_type='fast'
  character(len=20) :: comms_arch='serial'
  character(len=20) :: source_time=''
  character(len=20) :: source_date=''
  character(len=20) :: compile_date='Thu 16 Jul 2020'
  character(len=20) :: compile_time='11:48 CEST'
 end type build_info_type
 type(build_info_type), public, save :: build_info
endmodule od_build
