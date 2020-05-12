BEGIN{
    count=0
    OFS=","
    print id""count,id"_name","sequence",id"_len" 
    # out for fasta
    printf "" > out
    close(out)
}
/^>/{
    count++;
    name=substr($1,2);
    getline;
    print id""count,name,$1,length($1)
    printf ">%s\n%s\n",id""count,$1 >> out
    close(out)
}
