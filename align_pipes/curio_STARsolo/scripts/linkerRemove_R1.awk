(FNR-1) % 2 == 0 { name=$1; chr=$2; len=$3; next }
(FNR-2) % 4 == 0 { seq=substr($0,s,S) }
                 { print name "." seq, chr, len
                   print substr($0,n+1) }
