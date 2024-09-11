awk '/^>/ { # header pattern detected
        if (seqlen){
         # print previous seqlen if exists
         print seqlen
         }

         # print the tag
         print

         # initialize sequence
         seqlen = 0

         # skip further processing
         next
      }

# accumulate sequence length
{
seqlen += length($0)
}
# remnant seqlen if exists
END{if(seqlen){print seqlen}}' \
| paste -d " " - - \
| tr -d ">" \
| awk 'BEGIN{OFS="\t"}{ print $1,$2 }'
