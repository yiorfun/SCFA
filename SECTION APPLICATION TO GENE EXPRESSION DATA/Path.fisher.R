Path.fisher <- function(x, background, pathway, top_num) {
  ####x is the list of query genes
  ####backgroud is a list of background genes that query genes from 
  ####pathway is a list of different pathway genes
  
  count_table<-matrix(0,2,2)
  
  x<-toupper(x)
  background<-toupper(background)
  
  index<-which(background %in% x==FALSE)
  background_non_gene_list<-background[index]
  
  pathway<-lapply(pathway,function(x) intersect(background,x))
  
  get.fisher <- function(path) {
    res <- NA

    ####in the gene list and in the pathway
    count_table[1,1]<-sum(x %in% path)

    ####in the gene list but not in the pathway
    count_table[1,2]<-length(x)-count_table[1,1]

    ####not in the gene list but in the pathway
    count_table[2,1]<-sum(background_non_gene_list%in% path)

    ####not in the gene list and not in the pathway
    count_table[2,2]<-length(background_non_gene_list)-count_table[2,1]       

    matched_gene<-x[x %in% path]
 
    match_num<-length(matched_gene)

    overlap_info<-array(0,dim=4)
    names(overlap_info)<-c("DE in Geneset","DE not in Genese","NonDE in Geneset","NonDE out of Geneset")
    overlap_info[1]=count_table[1,1]
    overlap_info[2]=count_table[1,2]
    overlap_info[3]=count_table[2,1]
    overlap_info[4]=count_table[2,2]
    if(length(count_table)==4){
      res <- fisher.test(count_table, alternative="greater")$p}
    return(list(p_value=res,match_gene=matched_gene,match_num=match_num,
                fisher_table=overlap_info))
  }
  p_val<-array(0,dim=length(pathway))
  
  match_gene_list<-list(length(pathway))
  num1<-array(0,dim=length(pathway))
  num2<-matrix(0,nrow=length(pathway),ncol=4)
  colnames(num2)<-c("DE in Geneset","DE not in Genese","NonDE in Geneset","NonDE out of Geneset")
  
  for(i in 1:length(pathway)){
    result<-get.fisher(pathway[[i]])
    p_val[i]<-result$p_value
    match_gene_list[[i]]<-result$match_gene
    num1[i]<-result$match_num
    num2[i,]<-result$fisher_table
  }
  names(p_val) <-names(pathway)
  q_val <- p.adjust(p_val, "BH")
  
  sig_index<-sort(p_val,decreasing=FALSE,index.return=TRUE,na.last=NA)$ix  

  sig_index_top<-sig_index[1:top_num]
  match_gene_fraction<-array(0,dim=top_num)
  match_gene<-array(0,dim=top_num)
  for(i in 1:top_num){  
    match_gene[i]<-paste(match_gene_list[[sig_index_top[i]]],collapse="/")    
  }   
  
  summary<-format(data.frame(pvalue=p_val[sig_index_top],
                      qvalue=q_val[sig_index_top],
                      DE_in_Set=num2[sig_index_top,1],
                      DE_not_in_Set=num2[sig_index_top,2],
                      NonDE_in_Set=num2[sig_index_top,3],
                      NonDE_not_in_Set=num2[sig_index_top,4],
                      Odds_Ratio=(num2[sig_index_top,1]*num2[sig_index_top,4])/(num2[sig_index_top,2]*num2[sig_index_top,3]),
                      log_Odds_Ratio=log((num2[sig_index_top,1]*num2[sig_index_top,4])/(num2[sig_index_top,2]*num2[sig_index_top,3])),
                      Pathway_Size=num2[sig_index_top,1]+num2[sig_index_top,3],
                      match_gene=match_gene),
                      digits=3)
                          
  
  return(summary)
}

