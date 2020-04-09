save.table=function(x, ofn, verbose=T)
{
    if (verbose) cat(sprintf("saving result table: %s\n", ofn))
    write.table(x, ofn, quote=F, row.names=F, sep="\t")
}

# seq.vec: list of strings
# xcoverage: vector of read xcoverage (same length as seq.vec)
# ofn: filename
# line.width: length of line

save.seq.vector=function(seq.list, xcov, circ, ofn.fasta, ofn.table, line.width=80)
{
    df = data.frame(contig=names(seq.list), length=0, xcov=xcov, circ=circ)
    N = length(seq.list)
    lines = NULL
    for (i in 1:N) {
        contig = names(seq.list)[i]
        seq.v = seq.list[[i]]
        M = length(seq.v)
        mm = t(matrix("", ceiling(M/line.width), line.width))
        mm[1:M] = seq.v
        mm = t(mm)
        ss = apply(mm, 1, function(x) paste(x, collapse=""))
        lines = c(lines, sprintf(">%s", contig))
        lines = c(lines, ss)
        df$length[i] = M
    }

    cat(sprintf("saving fasta file: %s\n", ofn.fasta))
    fc = file(ofn.fasta)
    writeLines(lines, fc)
    close(fc)

    save.table(df, ofn.table)
}

write.ref.fa=function(names,seqs,ofn.fasta) {
    cat("writing reference fasta...\n")
    dir = unlist(strsplit(ofn.fasta,split='/'))
    dir = paste(dir[1:length(dir)-1],collapse="/")
    ofn.fa = paste(dir,"ref.fa",sep='/')
    num.seqs = length(names)
    cat(length(seqs))
    lines = NULL
    for (i in 1:num.seqs){
        lines = c(lines,sprintf(">%s",names[[i]]))
	lines = c(lines,seqs[[i]])
    }
    f = file(ofn.fa)
    writeLines(lines,f)
    close(f)
}

get.seq=function(N) {
    nts = c("A", "G", "C", "T")
    seq.i = floor(runif(N) * 4) + 1
    seq.v = nts[seq.i]
}

get.ind=function(N) {
    seq.i = floor(runif(N)*4) + 1
}

get.forw=function(seq.i){
    nts = c("A", "G", "C", "T")
    seq.v = nts[seq.i]
}

get.com_rev=function(seq.i){
    nts = c("T","C","G","A")
    seq.v = nts[rev(seq.i)]
}

get.mut=function(seq.i,N){
    seq.mut = seq.i
    mut_coords = sample.int(length(seq.i),N)
    seq.mut[mut_coords] = floor(runif(N)*4) + 1
    seq.mut
}

gen.three_rep=function(ofn.fasta, ofn.table){
  plas1 = get.seq(20000)
  rep = get.seq(1000)
  plas2 = get.seq(2000)
  seq.list = list(m1=c(plas1,rep,rep),m2=c(rep))
  xcov = c(100,100)
  circ = c(T,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}


gen.three_tandem_rep=function(ofn.fasta, ofn.table){
  plas1 = get.seq(20000)
  rep = get.seq(1000)
  seq.list = list(m1=c(plas1,rep,rep,rep))
  xcov = c(100)
  circ = c(T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.10mutations=function(ofn.fasta, ofn.table)
{
  plas1 = get.seq(20000)
  rep_index = get.ind(1000)
  rep_1 = get.forw(rep_index)
  list_rep = NULL
  for (i in 1:10){
    list_rep = append(list_rep,list(get.forw(get.mut(rep_index,100))))
  }
  seq.list = list(m1=c(plas1,rep_1), m2=c(list_rep[[1]]),m3=c(list_rep[[2]]),m4=c(list_rep[[3]]),m5=c(list_rep[[4]]),m6=c(list_rep[[5]]),m7=c(list_rep[[6]]),m8=c(list_rep[[7]]),m9=c(list_rep[[8]]),m10=c(list_rep[[9]]))
  xcov = c(200,10,10,10,10,10,10,10,10,10)
  circ = c(T,F,F,F,F,F,F,F,F,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.10_mutations_x10=function(ofn.fasta, ofn.table)
{
  plas1 = get.seq(20000)
  rep_index = get.ind(1000)
  rep_1 = get.forw(rep_index)
  list_rep = vector("list",10)
  for (i in 1:10){
    list_rep[[i]] = list(get.forw(get.mut(rep_index,10)))
  }
  seq.list = vector("list",10)
  names = vector("list",10)
  seq.list[[1]] = c(plas1,rep_1)
  names[[1]] = 'm1'
  for (i in 2:10){
    names[[i]] = paste("m",toString(i),sep='')
    seq.list[[i]] = list(c(list_rep[[i-1]]))
  }
  names(seq.list) = names
  xcov = c(200,10,10,10,10,10,10,10,10,10)
  circ = c(T,F,F,F,F,F,F,F,F,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.100_mutations_x100=function(ofn.fasta, ofn.table)
{
  plas1 = get.seq(20000)
  rep_index = get.ind(1000)
  rep_1 = get.forw(rep_index)
  list_rep = vector("list",100)
  seq.list = vector("list",100)
  names = vector("list",100)
  for (i in 1:100){
    list_rep[[i]] = list(get.forw(get.mut(rep_index,100)))
  }
  seq.list[[1]] = c(plas1,rep_1)
  names[[1]] = 'm1'
  for (i in 2:100){
    names[[i]] = paste("m",toString(i),sep='')
    seq.list[[i]] = list(c(list_rep[[i-1]]))
  }
  names(seq.list) = names
  xcov = c(200,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10)
  circ = c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.1000_mutations_x100=function(ofn.fasta, ofn.table)
{
  plas1 = get.seq(20000)
  rep_index = get.ind(1000)
  rep_1 = get.forw(rep_index)
  list_rep = vector("list",100)
  seq.list = vector("list",100)
  names = vector("list",100)
  for (i in 1:100){
    list_rep[[i]] = list(get.forw(get.mut(rep_index,1000)))
  }
  seq.list[[1]] = c(plas1,rep_1)
  names[[1]] = 'm1'
  for (i in 2:100){
    names[[i]] = paste("m",toString(i),sep='')
    seq.list[[i]] = list(c(list_rep[[i]]))
  }
  names(seq.list) = names
  print(seq.list[[3]])
  print(seq.list[[2]])
  print(names)
  xcov = c(200)
  circ = c(T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.200_mutations_x100=function(ofn.fasta, ofn.table)
{
  plas1 = get.seq(20000)
  rep_index = get.ind(1000)
  rep_1 = get.forw(rep_index)
  list_rep = NULL
  for (i in 1:100){
    list_rep = append(list_rep,(get.forw(get.mut(rep_index,200))))
  }
  seq.list = NULL
  names = NULL
  seq.list = append(seq.list, list(c(plas1,rep_1)))
  names  = append(names,'m1')
  for (i in 2:100){
    names = append(names,paste("m",toString(i),sep=''))
    rep = list_rep[((i-1)*1000+1):(i*1000)]
    seq.list = append(seq.list, list(rep))
  }
  names(seq.list) = names
  xcov = c(200,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10)
  circ = c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}


gen.10_mutations_x100=function(ofn.fasta, ofn.table)
{
  plas1 = get.seq(20000)
  rep_index = get.ind(1000)
  rep_1 = get.forw(rep_index)
  list_rep = NULL
  for (i in 1:100){
    list_rep = append(list_rep,(get.forw(get.mut(rep_index,10))))
  }
  seq.list = NULL
  names = NULL
  seq.list = append(seq.list, list(c(plas1,rep_1)))
  names  = append(names,'m1')
  for (i in 2:100){
    names = append(names,paste("m",toString(i),sep=''))
    rep = list_rep[((i-1)*1000+1):(i*1000)]
    seq.list = append(seq.list, list(rep))
  }
  names(seq.list) = names
  xcov = c(200,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10)
  circ = c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}


# gen.10_mutations_x10=function(ofn.fasta, ofn.table)
# {
#   plas1 = get.seq(20000)
#   rep_index = get.ind(1000)
#   rep_1 = get.forw(rep_index)
#   list_rep = NULL
#   for (i in 1:10){
#     list_rep = append(list_rep,list(get.forw(get.mut(rep_index,10))))
#   }
#   for (i in 1:100){
#     append(seq)
#   }
#   seq.list
#   seq.list = list(m1=c(plas1,rep_1), m2=c(list_rep[[1]]),m3=c(list_rep[[2]]),m4=c(list_rep[[3]]),m5=c(list_rep[[4]]),m6=c(list_rep[[5]]),m7=c(list_rep[[6]]),m8=c(list_rep[[7]]),m9=c(list_rep[[8]]),m10=c(list_rep[[9]]))
#   xcov = c(200,10,10,10,10,10,10,10,10,10)
#   circ = c(T,F,F,F,F,F,F,F,F,F)
#   save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
# }

gen.70_mutations_x10=function(ofn.fasta, ofn.table)
{
  plas1 = get.seq(20000)
  rep_index = get.ind(1000)
  rep_1 = get.forw(rep_index)
  list_rep = NULL
  for (i in 1:10){
    list_rep = append(list_rep,(get.forw(get.mut(rep_index,70))))
  }
  seq.list = NULL
  names = NULL
  seq.list = append(seq.list, list(c(plas1,rep_1)))
  names  = append(names,'m1')
  for (i in 2:10){
    names = append(names,paste("m",toString(i),sep=''))
    rep = list_rep[((i-1)*1000+1):(i*1000)]
    seq.list = append(seq.list, list(rep))
  }
  names(seq.list) = names
  xcov = c(200,10,10,10,10,10,10,10,10,10)
  circ = c(T,F,F,F,F,F,F,F,F,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.400_mutations_x10=function(ofn.fasta, ofn.table)
{
  plas1 = get.seq(20000)
  rep_index = get.ind(1000)
  rep_1 = get.forw(rep_index)
  list_rep = NULL
  for (i in 1:10){
    list_rep = append(list_rep,(get.forw(get.mut(rep_index,400))))
  }
  seq.list = NULL
  names = NULL
  seq.list = append(seq.list, list(c(plas1,rep_1)))
  names  = append(names,'m1')
  for (i in 2:10){
    names = append(names,paste("m",toString(i),sep=''))
    rep = list_rep[((i-1)*1000+1):(i*1000)]
    seq.list = append(seq.list, list(rep))
  }
  names(seq.list) = names
  xcov = c(200,10,10,10,10,10,10,10,10,10)
  circ = c(T,F,F,F,F,F,F,F,F,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}


gen.tandem_same_rep=function(ofn.fasta, ofn.table)
{
  a = get.seq(10000)
  b = get.seq(1000)
  seq.list = list(m1=c(a,b,b))
  xcov = c(100)
  circ = c(T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}



# gen.70_mutations_x10=function(ofn.fasta, ofn.table)
# {
#   plas1 = get.seq(20000)
#   rep_index = get.ind(1000)
#   rep_1 = get.forw(rep_index)
#   list_rep = NULL
#   for (i in 1:10){
#     list_rep = append(list_rep,list(get.forw(get.mut(rep_index,70))))
#   }
#   seq.list = list(m1=c(plas1,rep_1), m2=c(list_rep[[1]]),m3=c(list_rep[[2]]),m4=c(list_rep[[3]]),m5=c(list_rep[[4]]),m6=c(list_rep[[5]]),m7=c(list_rep[[6]]),m8=c(list_rep[[7]]),m9=c(list_rep[[8]]),m10=c(list_rep[[9]]))
#   typeof(seq.list)
#   xcov = c(200,10,10,10,10,10,10,10,10,10)
#   circ = c(T,F,F,F,F,F,F,F,F,F)
#   save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
# }




gen.30_mutations_x10=function(ofn.fasta, ofn.table)
{
  plas1 = get.seq(20000)
  rep_index = get.ind(1000)
  rep_1 = get.forw(rep_index)
  list_rep = NULL
  for (i in 1:10){
    list_rep = append(list_rep,list(get.forw(get.mut(rep_index,30))))
  }
  seq.list = list(m1=c(plas1,rep_1), m2=c(list_rep[[1]]),m3=c(list_rep[[2]]),m4=c(list_rep[[3]]),m5=c(list_rep[[4]]),m6=c(list_rep[[5]]),m7=c(list_rep[[6]]),m8=c(list_rep[[7]]),m9=c(list_rep[[8]]),m10=c(list_rep[[9]]))
  xcov = c(200,10,10,10,10,10,10,10,10,10)
  circ = c(T,F,F,F,F,F,F,F,F,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.reverse=function(ofn.fasta, ofn.table)
{
  plas1 = get.seq(20000)
  plas2 = get.seq(20000)
  rep_index = get.ind(1000)
  a = get.seq(5000)
  b = get.seq(5000)
  rep_forw = get.forw(rep_index)
  rep_rev = get.com_rev(rep_index)
  seq.list = list(m1=c(plas1,rep_forw,a), m2=c(plas2,rep_rev,b))
  names(seq.list)
  xcov = c(200,400)
  circ = c(T,T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.tan_rep_rev_2plas=function(ofn.fasta, ofn.table)
{
  plas1 = get.seq(20000)
  plas2 = get.seq(20000)
  rep_index = get.ind(1000)
  rep_forw = get.forw(rep_index)
  rep_rev = get.com_rev(rep_index)
  rep_index2 = get.ind(1000)
  rep_forw2 = get.forw(rep_index2)
  rep_rev2 = get.com_rev(rep_index2)
  seq.list = list(m1=c(plas1,rep_forw,rep_forw2), m2=c(plas2,rep_rev,rep_rev2))
  xcov = c(200,50)
  circ = c(T,T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.tandemreverserep=function(ofn.fasta, ofn.table)
{
  plas1 = get.seq(20000)
  plas2 = get.seq(20000)
  rep_index = get.ind(1000)
  rep_forw = get.forw(rep_index)
  rep_rev = get.com_rev(rep_index)
  seq.list = list(m1=c(plas1,rep_forw,rep_rev), m2=c(plas2,rep_forw))
  xcov = c(200,50)
  circ = c(T,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.2repreverse=function(ofn.fasta, ofn.table)
{
  plas1 = get.seq(20000)
  plas2 = get.seq(20000)
  rep_index = get.ind(1000)
  a = get.seq(5000)
  b = get.seq(5000)
  rep_forw = get.forw(rep_index)
  rep_rev = get.com_rev(rep_index)
  rep2_index = get.ind(1000)
  rep2_forw = get.forw(rep2_index)
  rep2_rev = get.com_rev(rep2_index)
  seq.list = list(m1=c(plas1,rep_forw,a,rep2_rev), m2=c(plas2,rep_rev,b,rep2_forw))
  xcov = c(200,400)
  circ = c(T,T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.example=function(ofn.fasta, ofn.table)
{
    a = get.seq(1000)
    b = get.seq(10000)
    p = get.seq(2000)
    seq.list = list(m1=c(a,p,b), m2=p)
    xcov = c(100,1000)
    circ = c(F,T)
    save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.lin=function(ofn.fasta, ofn.table)
{
    a = get.seq(100000)
    b = get.seq(100)
    seq.list = list(contig1 = a, contig2 = b)
    xcov = c(100,100)
    circ = c(F,F)
    save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.plas1=function(ofn.fasta, ofn.table)
{
    plas1=get.seq(50000)
    seq.list=list(m1= plas1)
    xcov=1000
    circ=T
    save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.2plas=function(ofn.fasta, ofn.table)
{
  plas1=get.seq(50000)
  plas2=get.seq(25000)
  seq.list=list(m1=plas1,me=plas2)
  xcov=c(400,200)
  circ=c(T,T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}


gen.repPlasChrom=function(ofn.fasta, ofn.table)
{
    rep = get.seq(600)
    a = get.seq(10000)
    b = get.seq(17000)
    plas1a = get.seq(10000)
    plas1b = get.seq(10000)
    seq.list = list(m1=c(a,rep,b), m2=c(plas1a,rep,plas1b))
    xcov = c(1000, 1000)
    circ = c(F,T)
    save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.rep2Plas=function(ofn.fasta, ofn.table)
{
  rep = get.seq(1000)
  plas1 = get.seq(10000)
  plas2 = get.seq(20000)
  seq.list = list(m1=c(plas1,rep), m2=c(plas2,rep))
  xcov = c(100, 250)
  circ = c(T, T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.rep2Plasedit=function(ofn.fasta, ofn.table)
{
  rep = get.seq(1147)
  plas1 = get.seq(10000)
  plas2 = get.seq(20000)
  seq.list = list(m1=c(plas1,rep), m2=c(plas2,rep))
  xcov = c(300, 250)
  circ = c(T, T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.2repPlasChrom=function(ofn.fasta, ofn.table)
{
  rep1 = get.seq(600)
  rep2 = get.seq(1000)
  a = get.seq(10000)
  b = get.seq(17000)
  c = get.seq(10000)
  plas1a = get.seq(10000)
  plas1b = get.seq(20000)
  plas1c = get.seq(15500)
  seq.list = list(m1=c(a,rep1,b,rep2,c), m2=c(plas1a,rep1,plas1c))
  xcov = c(1000, 1000)
  circ = c(F,T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.rep3PlasChrom=function(ofn.fasta, ofn.table)
{
  rep = get.seq(800)
  plas1a = get.seq(10000)
  plas1b = get.seq(15000)
  plas2a = get.seq(9000)
  plas2b = get.seq(18000)
  a = get.seq(10000)
  b = get.seq(9000)
  seq.list = list(m1=c(plas1a,rep,plas1b), m2=c(plas2a,rep,plas2b), m3=c(a,rep,b))
  xcov = c(100,200,40)
  circ = c(T,T,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.2plas2chrom=function(ofn.fasta, ofn.table)
{
  plas1 = get.seq(10000)
  plas2 = get.seq(10000)
  chrom1 = get.seq(10000)
  chrom2 = get.seq(10000)
  seq.list = list(m1=c(plas1), m2=c(plas2), m3=c(chrom1), m4=c(chrom2))
  xcov = c(100,100,100,100)
  circ = c(T,T,F,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.rep2PlasChrom=function(ofn.fasta, ofn.table)
{
  rep = get.seq(600)
  plas1a = get.seq(5000)
  plas1b = get.seq(20000)
  plas2a = get.seq(7000)
  plas2b = get.seq(25000)
  chrom = get.seq(10000)
  seq.list = list(m1=c(plas1a,rep,plas1b), m2=c(plas2a,rep,plas2b), m3=c(chrom))
  xcov = c(1000, 900, 400)
  circ = c(T, T, F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.2rep2plas=function(ofn.fasta, ofn.table)
{
  rep1 = get.seq(1000)
  rep2 = get.seq(1000)
  plas1a = get.seq(10000)
  plas1b = get.seq(17000)
  plas1c = get.seq(10000)
  plas2a = get.seq(10000)
  plas2b = get.seq(20000)
  plas2c = get.seq(15500)
  seq.list = list(m1=c(plas1a,rep1,plas1b,rep2,plas1c), m2=c(plas2a,rep1,plas2b,rep2,plas2c))
  xcov = c(100, 100)
  circ = c(T,T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.rep2Plasrepchrom=function(ofn.fasta, ofn.table)
{
  rep1 = get.seq(600)
  rep2 = get.seq(1000)
  plas1a = get.seq(10000)
  plas1b = get.seq(17000)
  plas2a = get.seq(10000)
  plas2b = get.seq(20000)
  plas2c = get.seq(15500)
  a = get.seq(7000)
  seq.list = list(m1=c(plas1a,rep1,plas1b), m2=c(plas2a,rep1,plas2b,rep2,plas2c), m3=c(a,rep2))
  xcov = c(1000, 1000, 1000)
  circ = c(T,T,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.rep3Plas=function(ofn.fasta, ofn.table)
{
  rep = get.seq(600)
  plas1 = get.seq(10000)
  plas2 = get.seq(10000)
  plas3 = get.seq(12000)
  seq.list = list(m1=c(plas1,rep), m2=c(plas2,rep), m3=c(plas3,rep))
  xcov = c(1000, 1000, 1000)
  circ = c(T,T,T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.2Plas=function(ofn.fasta, ofn.table)
{
  plas1 = get.seq(10000)
  plas2 = get.seq(17000)
  seq.list = list(m1=c(plas1), m2=c(plas2))
  xcov = c(1000, 1000)
  circ = c(T,T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.repPlas3chrom=function(ofn.fasta, ofn.table)
{
  plas = get.seq(10000)
  rep = get.seq(1000)
  chrom1 = get.seq(15000)
  chrom2 = get.seq(8000)
  chrom3 = get.seq(12000)
  seq.list = list(m1=c(plas,rep), m2=c(chrom1,rep), m3=c(chrom2,rep), m4=c(chrom3,rep))
  xcov = c(1000, 1000, 1000, 1000)
  circ = c(T,F,F,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.2rep4plas=function(ofn.fasta, ofn.table)
{
  rep1 = get.seq(900)
  rep2 = get.seq(600)
  plas1 = get.seq(10000)
  plas2 = get.seq(15000)
  plas3 = get.seq(12000)
  plas4 = get.seq(11000)
  seq.list = list(m1=c(plas1,rep1), m2=c(plas2,rep1), m3=c(plas3,rep2), m4=c(plas4,rep2))
  xcov = c(1000, 1000, 1000, 1000)
  circ = c(T,T,T,T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.plaschrom=function(ofn.fasta, ofn.table)
{
  plas = get.seq(10000)
  chrom = get.seq(50000)
  seq.list = list(m1=c(plas), m2=c(chrom))
  xcov = c(1000, 1000)
  circ = c(T,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.3plas2rep=function(ofn.fasta, ofn.table){
  rep1 = get.seq(900)
  rep2 = get.seq(600)
  plas1 = get.seq(10000)
  plas2 = get.seq(15000)
  plas2b = get.seq(1000)
  plas3 = get.seq(12000)
  seq.list = list(m1=c(plas1,rep1), m2=c(plas2,rep1,plas2b,rep2), m3=c(plas3,rep2))
  xcov = c(1000, 1000, 1000)
  circ = c(T,T,T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.rep_plas_chrom=function(ofn.fasta, ofn.table){
  rep = get.seq(10000)
  plas1 = get.seq(10000)
  plas2 = get.seq(10000)
  seq.list = list(m1=c(plas1,rep), m2=c(plas2,rep))
  xcov = c(200, 200)
  circ = c(T,T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.6plas=function(ofn.fasta, ofn.table){
  plas1 = get.seq(1500)
  plas2 = get.seq(5000)
  plas3 = get.seq(10000)
  plas4 = get.seq(20000)
  plas5 = get.seq(40000)
  plas6 = get.seq(75000)
  seq.list = list(m1=c(plas1), m2=c(plas2), m3=c(plas3), m4=c(plas4), m5=c(plas5), m6=c(plas6))
  xcov = c(100,100,100,100,100,100)
  circ = c(T,T,T,T,T,T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.dilute6plas=function(ofn.fasta, ofn.table){
  plas1 = get.seq(15000)
  plas2 = get.seq(15000)
  plas3 = get.seq(15000)
  plas4 = get.seq(15000)
  plas5 = get.seq(15000)
  plas6 = get.seq(15000)
  seq.list = list(m1=c(plas1), m2=c(plas2), m3=c(plas3), m4=c(plas4), m5=c(plas5), m6=c(plas6))
  xcov = c(500,250,100,50,25,10)
  circ = c(T,T,T,T,T,T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.repplaschrom=function(ofn.fasta, ofn.table){
  plas1 = get.seq(15000)
  rep = get.seq(2000)
  chrom = get.seq(1000)
  seq.list = list(m1=c(plas1,rep), m2=c(rep,chrom))
  xcov = c(100,500)
  circ = c(T,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.chrom_rep=function(ofn.fasta, ofn.table) {
  rep = get.seq(500)
  chrom1 = get.seq(10000)
  chrom2 = get.seq(10000)
  chrom3 = get.seq(1000)
  seq.list = list(m1=c(chrom1,rep,chrom2,rep,chrom3))
  xcov = c(50)
  circ = c(F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.linear_tandem=function(ofn.fasta, ofn.table) {
  rep = get.seq(1000)
  seq.list = list(m1=c(rep,rep,rep,rep))
  xcov = c(20)
  circ = c(F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.type1plas2=function(ofn.fasta, ofn.table) {
  plas = get.seq(15000)
  rep = get.seq(5000)
  chrom1 = get.seq(10000)
  chrom2 = get.seq(10000)
  seq.list = list(m1=c(plas,rep), m2=c(chrom1,rep,chrom2))
  xcov = c(200,1000)
  circ = c(T,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.sim_split=function(ofn.fasta, ofn.table) {
  rep = get.seq(5000)
  a = get.seq(10000)
  b = get.seq(10000)
  seq.list = list(m1=c(rep,a), m2=c(rep,b))
  xcov = c(200,400)
  circ = c(F,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.p1=function(ofn.fasta, ofn.table) {
  a = get.seq(10000)
  seq.list = list(m1=c(a))
  xcov = c(100)
  circ = c(T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.type1plas=function(ofn.fasta, ofn.table) {
  plas = get.seq(15000)
  rep = get.seq(5000)
  chrom1 = get.seq(10000)
  chrom2 = get.seq(10000)
  seq.list = list(m1=c(plas,rep), m2=c(chrom1,rep,chrom2))
  write.ref.fa(c("plas","rep","chrom1","chrom2"),c(paste(plas,collapse=""),paste(rep,collapse=""),paste(chrom1,collapse=""),paste(chrom2,collapse="")),ofn.fasta)
  xcov = c(100,500)
  circ = c(T,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.type2plas=function(ofn.fasta, ofn.table) {
  plas = get.seq(15000)
  rep = get.seq(5000)
  chrom1 = get.seq(10000)
  chrom2 = get.seq(10000)
  seq.list = list(m1=c(plas,rep), m2=c(chrom1,rep,chrom2))
  write.ref.fa(c("plas","rep","chrom1","chrom2"),c(paste(plas,collapse=""),paste(rep,collapse=""),paste(chrom1,collapse=""),paste(chrom2,collapse="")),ofn.fasta)
  xcov = c(500,100)
  circ = c(T,F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.tandem=function(ofn.fasta, ofn.table) {
  rep = get.seq(5000)
  a = get.seq(10000)
  b = get.seq(10000)
  seq.list = list(m1=c(a,rep,rep,b))
  write.ref.fa(c("a","rep","b"),c(paste(a,collapse=""),paste(rep,collapse=""),paste(b,collapse="")),ofn.fasta)
  xcov = c(200)
  circ = c(F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.figure_eight=function(ofn.fasta, ofn.table) {
  plas1 = get.seq(20000)
  rep = get.seq(5000)
  plas2 = get.seq(10000)
  seq.list = list(m1=c(plas1,rep), m2=c(plas2,rep))
  write.ref.fa(c("plas1","rep","plas2"),c(paste(plas1,collapse=""),paste(rep,collapse=""),paste(plas2,collapse="")),ofn.fasta)
  xcov = c(100,50)
  circ = c(T,T)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

gen.rep_chrom=function(ofn.fasta, ofn.table) {
  rep = get.seq(5000)
  a = get.seq(10000)
  b = get.seq(10000)
  c = get.seq(10000)
  seq.list = list(m1=c(a,rep,b,rep,c))
  write.ref.fa(c("a","rep","b","c"),c(paste(a,collapse=""),paste(rep,collapse=""),paste(b,collapse=""),paste(c,collapse="")),ofn.fasta)
  xcov = c(100)
  circ = c(F)
  save.seq.vector(seq.list=seq.list, xcov=xcov, circ=circ, ofn.fasta=ofn.fasta, ofn.table=ofn.table)
}

