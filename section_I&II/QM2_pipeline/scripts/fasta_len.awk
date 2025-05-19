#!/usr/bin/awk -f

BEGIN {
   OFS = "\t";
   #!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }
}

{ 
  # first sequence
  if ( $0 ~ "^>" && FNR == 1 ) {
      id = substr($0,2)
  # other sequences
  } else if ($0 ~ "^>" && FNR != 1) {
      # prints previous sequence lengh and resets len var
      print id, len
      len = 0
      # stores new id
      id = substr($0,2)
  # nucleotide sequence
  } else if ( $0 !~ "^>" ) {
      len += length($0)
  }
}

# prints length of last sequence 
END { print id, len }
