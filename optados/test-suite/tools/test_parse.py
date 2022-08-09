# Quick wrapper to test a parser on a partiuclar file
# AJM 8/2/2022
#
import sys
import parsers.parse_od2od_fmt

print(sys.argv)
out=parsers.parse_od2od_fmt.parse(sys.argv[1])
print(out)
