#!/usr/bin/awk -f
BEGIN{
    OFS=","
    print "sequence","read_count"

}

/^>/{
    split($1,array,"_")
    getline;
    print $1,array[2];
}
