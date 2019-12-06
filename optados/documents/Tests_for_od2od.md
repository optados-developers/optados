#Test the Tester
AJM 6/12/19

Where RB, WB, RF, WF stand for read binary, write binary, read formatted
  and write formatted respectively.

# Test the Tester
We can assume that RB is in working order as it uses the same routines as
  OptaDOS. Hence the internal state of od2od can be assumed "correct" if
  od2od reads in a binary file, call it b1.

##Test 1
 b1 --> bin fmt --> f1
 b1 --> bin bin --> b2
 b2 --> bin fmt --> f2

Now if f1=f2, then since b1 and b2 had the same operation (bin fmt) applied
  to them, b1=b2.  od2od can faithfully remake a binary file.

##Test 2
f1 --> fmt fmt --> f3

Now if f1=f3 we can say that WF at least is the reciprocal of RF.  However 
  this does not prove that RF creates the correct internal state of od2od.
  However, since RB /does/ failthfully create the correct internal state,
  any formatted file made by od2od will faithfully create the correct 
  internal state when read back in. Since any error introduced by WF is 
  undone by RF even if we cannot prove that the formatted file is correct.

#Test suite.
This will store a formatted version of each of the files.

##TS1
f1 --> fmt fmt --> fs2
Then check that f1=f2 tests the RF, WF code paths.

##TS2
f1 --> fmt bin --> b1
b1 --> bin fmt --> f2
Then check that f1=f2 tests the RB, WB code paths.
