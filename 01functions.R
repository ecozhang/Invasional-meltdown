#----------------- for both plant and soil analyses -------------------#


logit <- function(n){
  t <- log(n / (1 - n))
  return(t)
}



# add sign of significant level and format of the table
lrt.table <- function(table){
  table$" " = ifelse(table$`Pr(>Chi)` < 0.05, '*',
                     ifelse(table$`Pr(>Chi)` < 0.1, '\u2020', ''))#add significant levels
  return(table %>% as.data.frame() %>% 
           kable(digits = 3))
}



# extract random effects
extract_rand <- function(model){
  t_rand <- VarCorr(model)[,2] %>% 
    as.data.frame() %>% 
    rename('SD' = ".") %>% 
    unique() %>%
    rownames_to_column('Random effects') %>% 
    mutate(    SD           = sprintf('%.3f',round(as.numeric(as.character(SD)),3))) %>% 
    mutate(    SD           = as.character(SD),
               `Random effects` = str_replace(`Random effects`, 'alone.*|[A-Z].*', ''))
  
  t_rand[nrow(t_rand),1] <- 'Residual'
  return(t_rand)
}



# extract random effects for nested random effects (i.e. pdblock not used)
# note, this is only for the soil analyses
extract_rand_nest <- function(model){
  t_rand <- VarCorr(model)[,2] %>% 
    as.data.frame() %>% 
    rename('SD' = ".") %>% 
    unique() %>%
    rownames_to_column('Random effects') %>% 
    mutate(    SD           = sprintf('%.3f',round(as.numeric(as.character(SD)),3))) %>% 
    mutate(    SD           = as.character(SD),
               `Random effects` = str_replace(`Random effects`, 'alone.*|[A-Z].*', '')) %>% 
    slice(2:4) %>% 
    mutate(`Random effects` = c('Family', 'Species', 'Residual'))
  
  return(t_rand)
}




# combine results of drop1 and random effects
merge_fix_rand <- function(table, model, nest = F){
  table_fix <- table %>% 
    as.data.frame() %>% 
    dplyr::select(-AIC, -Df) %>% # remove AIC and DF, which are not needed
    rownames_to_column('terms') %>% 
    mutate_if(is.numeric, function(x){sprintf('%.3f', round(x,3))}) %>%  # three digits
    column_to_rownames('terms')
  colnames(table_fix) <- c('X2', 'P') # make column names easier
  if (nest) {
    table_rand <- extract_rand_nest(model)
    } else {
      table_rand <- extract_rand(model)}
  
  length_fix <- nrow(table_fix)
  length_rand<- nrow(table_rand)
  table_fix[(length_fix+1):(length_fix + length_rand + 1),] <- ''
  rownames(table_fix)[(length_fix+1):(length_fix + length_rand +1)] <- c('Random effects',table_rand[,1]) # add terms of random effects
  table_fix[(length_fix+1),1] <- 'SD'
  table_fix[(length_fix+2):(length_fix + length_rand + 1),1] <- table_rand[,2]
  
  table_fix$P = ifelse(table_fix$P < 0.05   & table_fix$P != '', paste(table_fix$P, '*', sep = ''),
                       ifelse(table_fix$P < 0.1 & table_fix$P != '', paste(table_fix$P, '\u2020', sep = ''),
                              table_fix$P))#add significant levels
  return(table_fix)
}
#------------------------------------------------------------------#

#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#

#------------------------ only for plant analyses-------------------------------------#

# drop1 function: the models for plants that were grown alone

drop1_single <- function(model){
  model_r2 <- update(model,    .~. - (T1_empty + T2_home_away + T3_origin) : origin_p2) # main effects
  model_r3 <- update(model_r2, .~. -  T1_empty - T2_home_away - T3_origin -  origin_p2) ## covariate
  
  lambda_r1   <- drop1(model,     test="Chisq")[-c(1,2),] # two-way
  lambda_r2   <- drop1(model_r2,  test="Chisq")[-c(1,2),] # main
  lambda_r3   <- drop1(model_r3,  test="Chisq")[-1,]      # covariate
  
  table <- merge_fix_rand(rbind(lambda_r3, lambda_r2, lambda_r1), m_lambda)
  return(table)
}






# plot function, for Figure 3a and 3b
plot_alone <- function(data = dat_lambda, title = 'Biomass when grown alone', 
                       xlab = '\nSoil treatment', ylab = 'Aboveground biomass [g]\n', tag = 'a', 
                       variable = 'biomass', ylim = c(0, 3), legend = T
){
  at=c(1.2,1.8,# empty
       3.7,4.3,# own,i.e. same species
       6.2,6.8,# alien soil species
       8.7,9.3 # native soil species
  )
  colnames(data)[grepl(variable, colnames(data))] <- 'target_biomass'# change the name of target variable as 'biomass'.
  
  t_lambda <- data %>%
    # average of each sp_p1 x target combi
    dplyr::group_by(target, sp_p1, treat, origin_p2) %>%
    dplyr::summarise(mean_bio = mean(target_biomass)) %>%# unit as g
    # average of treat x origin of target
    group_by(treat, origin_p2)%>%
    dplyr::summarise(mean = mean(mean_bio), se = sd(mean_bio)/sqrt(n()))%>%
    arrange(match(treat, c("empty", "own", "alien", "native")))
  
  p <- t_lambda %>% ggplot() + 
    aes(x = at, y = mean, col = origin_p2) + 
    geom_point() + geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0) +
    scale_x_continuous(limits = c(0.5, 10), breaks = c(1.5,4,6.5,9), 
                       labels = c("non-conditioned", "home", "alien", "native")) + 
    scale_y_continuous(limits = ylim) +
    labs(x = xlab, y = ylab, title = title, col = 'Test species:', tag = tag) + 
    scale_color_manual(values = col_fig) +
    theme
  
  if (legend == F) p <- p + theme(legend.position = 'none') # add legend or not
  return(p)
}




#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#







#-----------------------------only for soil analyses ------------------------#

# calculate diversity, abundance for each its functinal group (type)
rel_abun<-function(otu_table,function_table,type){
  otu_table_sub<-otu_table[,as.character(function_table$OTU_ID[function_table[,type]%in%type])]
  abundance<-sum(otu_table_sub)/sum(otu_table)
  diversity<-ncol(otu_table_sub)/ncol(otu_table)
  print(paste('% abundance of',type,'=',abundance))
  print(paste('# diversity of',type,'=',ncol(otu_table_sub)))
  print(paste('% diversity of',type,'=',diversity))
}




# beta-diversity for each pair of species
beta_by_sp <- function(otu_table = otu_table_16s_nonster,
                       info = info_16s){
  
  beta_mat <- as.matrix(vegdist(otu_table)) 
  beta_mat[lower.tri(beta_mat)] <- 999 #lower matrix to 999, and later on remove it, as the matrix is symmetric
  diag(beta_mat) <- 999 # diagnol to 999, and later on remove it
  
  # note whether logit transform beta diversity
  table <- as.data.frame(beta_mat)%>%
    rownames_to_column(var = 'sample1') %>%
    gather(sample2, beta, -sample1) %>%
    filter(beta != 999)%>% # remove lower matrix and diagnol
    #--------add information of species 1 and species 2 ----------#
    left_join(info, by = c('sample1' ='sample_id'))%>%#add info for sp1
    dplyr::rename(status1 = status, ster.not1 = ster.not,species1 = species)%>%
    left_join(info, by = c('sample2' ='sample_id'))%>%#add info for sp2
    dplyr::rename(status2 = status, ster.not2 = ster.not,species2 = species)%>%
    mutate(species1 = as.character(species1),
           species2 = as.character(species2),
           status1 = as.character(status1),
           status2 = as.character(status2)) %>%
    mutate(combi = ifelse(species1 > species2, paste(species1, species2), paste(species2, species1)))%>%# get combination, and with order
    mutate(combi_status = ifelse(status1 > status2, paste(status1, status2, sep = '-'), paste(status2, status1, sep = '-'))) %>% 
    mutate(combi_status = ifelse(species1 == species2, paste('own', status1, sep = '-'), as.character(combi_status)))
  
  # group and average
  table <- table %>% 
    group_by(combi, combi_status)%>%
    summarise(mean_beta = mean(beta))%>% 
    filter((combi_status %in% c('native-native','native-alien','alien-alien', 'own-alien', 'own-native')))%>%
    separate(combi,sep = ' ',c('sp1','sp2')) %>% 
    # ------------add family------------#
    left_join(info_family %>% dplyr::select(species, family), by = c('sp1' = 'species')) %>% 
    rename(family1 = family) %>% 
    left_join(info_family %>% dplyr::select(species, family), by = c('sp2' = 'species')) %>% 
    rename(family2 = family) %>% 
    mutate(combi_status = as.factor(combi_status))
  
  # add contrasts
  table <- table %>% 
    mutate(intra_vs_inter = ifelse(combi_status %in% c('own-alien', 'own-native'), 0.6, -0.4),
           intra_a_vs_n = ifelse(combi_status == 'own-alien', 0.5,
                                 ifelse(combi_status == 'own-native', -0.5, 0)),
           aa_vs_na = ifelse(combi_status == 'alien-alien', 2/3,
                             ifelse(combi_status %in% c('native-alien', 'native-native'), -1/3, 0)),
           na_vs_nn = ifelse(combi_status == 'native-native', 2/3,
                             ifelse(combi_status %in% c('native-alien', 'alien-alien'), -1/3, 0)))


  return(table)
}





# get two matrices, one for the species-by-species phylogenetic matrix, 
# and another for the species-by-species soil-microbial dissimilarity
dist_mat <- function(beta_pd = beta_pd_its_nonster){
  number <- 10 # number of species
  mat_pd <- mat_beta <-matrix(,number,number)
  # arrange species according to status and species
  colnames(mat_pd) <-rownames(mat_pd) <- 
    colnames(mat_beta) <- rownames(mat_beta) <- (info_its %>%
                                                   filter(status %in% c('native','alien')) %>%
                                                   dplyr::select(species,status) %>%
                                                   distinct() %>%
                                                   arrange(status,species))[,'species']
  #generate two matrices, matrix beta and matrix pd, whose rows and columns represent species
  for (i in 1:number){
    row.i<-rownames(mat_beta)[i]
    for (j in 1:number){
      # for beta
      column.j <- colnames(mat_beta)[j]
      point  <-(beta_pd %>% filter(sp1==row.i,sp2==column.j))$mean_beta
      point2 <-(beta_pd %>% filter(sp1==column.j,sp2==row.i))$mean_beta
      value = c(point,point2)
      mat_beta[i,j]<-ifelse(length(value) ==0, NA, value)
      mat_beta[j,i]<-ifelse(length(value) ==0, NA, value)
      # for pd
      point  <-(beta_pd %>% filter(sp1==row.i,sp2==column.j))$pd
      point2 <-(beta_pd %>% filter(sp1==column.j,sp2==row.i))$pd
      value = c(point, point2)
      mat_pd[i,j]<-ifelse(length(value) ==0,NA,value)
      mat_pd[j,i]<-ifelse(length(value) ==0,NA,value)
    }
  }
  return(list(sqrt(mat_pd), logit(mat_beta))) # sqrt and logit transform 
}





#------------------------- plots for soil-----------------------#

# corrplot, figure 4a-d
plot_cor<-function(beta_mat_sp = beta_mat_sp_16s_nonster, 
                   title = 'Bacteria', 
                   number = 10, 
                   tag = 'a'){
  
  # get the sp x sp beta-diversity correlation matrix
  mat <- matrix(, number, number)
  colnames(mat) <- rownames(mat) <-(info_its %>% 
                                      filter(status %in% c('native','alien')) %>% 
                                      dplyr::select(species, status) %>%
                                      distinct() %>% 
                                      arrange(status,species))[,'species']
  # this loop is ugly, but it does its job. I could have made it more understandable.
  for (i in 1 : number){
    row.i <- rownames(mat)[i]
    for (j in 1:number){
      column.j <- colnames(mat)[j]
      point    <-(beta_mat_sp %>% filter(sp1 == row.i,    sp2 == column.j))$mean_beta
      point2   <-(beta_mat_sp %>% filter(sp1 == column.j, sp2 == row.i))$mean_beta
      value = c(point,point2)
      mat[i,j] <- ifelse(length(value) ==0, NA, value)
      mat[j,i] <- ifelse(length(value) ==0, NA, value)
    }
  }
  
  
  mat[lower.tri(mat)] <- NA # use half matrix, as it is symmetric
  # abbreviation of species name
  rownames(mat) <- colnames(mat)<- paste(str_extract(rownames(mat),'^\\w'), # first string of genus
        toupper(str_extract(str_extract(rownames(mat),'_\\w'),'[a-z]')), # first string of species and uppercase # space for lines
        sep = '')
  # center data
  mat2 <- (mat - mean(mat,na.rm=TRUE))
  mat2 <- mat2/max(abs(mat2),na.rm=TRUE)

  diag <- as.vector(diag(mat2))
  diag(mat2) <- NA
  mat4 <- rbind(diag, NA, mat2)
  mat4 <- cbind(NA, mat4)
  rownames(mat4)[1] <- ' '
  
  # plot
  corrplot(mat4,type="upper", diag = T, 
           na.label = " ", # how to show na
           mar = c(0.5, 2.5, 1, 0.1),
           method = 'color',
           col = colorRampPalette(c("red","white","blue"))(200),
           cl.pos = 'n', cl.ratio	= 0.13, cl.align = "r",cl.length = 5, # label
           tl.pos = 'n', tl.col= 'black',tl.cex = 0.5,#, tl.pos = 'n', tl.srt = 90,
           ) # tip
  
  #------------- additional information on the plot --------------#
  # title
  text(6.5, 14.5, title, font = 1, cex = 1)
  # species name
  text(2:11,    12.85,    colnames(mat4)[-1],   cex = 0.6, col = c(rep(col_fig[1],4), rep(col_fig[2],6)))
  text(2 + 0:8, 10 - 0:8, colnames(mat4)[2:10], cex = 0.6,  col = c(rep(col_fig[1],4), rep(col_fig[2],5))) # diagnol
  
  # border
  rect(1.5, 12.5, 11.5, 11.5, border  = 'gray25') # intra
  segments(5.5,12.5 ,5.5, 11.5, col = 'gray25') # seperate alien and native
  #inter
  segments((2:10+0.5), (10:2 + 0.5), (2:10 + 0.5), (10:2 - 0.5), lwd = 1.2)
  segments(2:10 + 0.5, 9:1 +0.5, 3:11 + 0.5, 9:1 + 0.5, lwd = 1.2)
  segments(2 + 0.5, 10 +0.5, 5 + 0.5, 10 + 0.5, lwd = 1.2)
  segments(11 + 0.5, 1 +0.5, 11 + 0.5, 6 + 0.5, lwd = 1.2)
  # native - alien
  rect(5.5, 6.5, 11.5, 10.5, lwd = 1.2)
  
  # tag
  mtext(tag, side = 3, line= -1.1, at= 0.5, cex = 0.7, font = 2)
}








# calculate abundance for each sample
abun_by_pot<-function(otu_table,info_table,total_abundance = 9500){
  table <- as.data.frame(apply(otu_table,1,sum))/total_abundance
  colnames(table)<-'rel_abu'
  table <- table %>%
    rownames_to_column(var = 'sample_id') %>%
    left_join(info_table)%>%
    mutate(fix = sub('.*ITS','',sample_id)) %>%
    left_join(info_family) %>% 
    #----- add contrast ----- #
    mutate(T1_empty = ifelse(status == 'non-conditioned', -2/3, 1/3),
           T2_alien_native = ifelse(status == 'non-conditioned', 0,
                                    ifelse(status == 'alien', 0.5, -0.5)))
  return(table)   
}

# calculate alpha diversity for each sample
alpha_div<-function(data,sample_info){
  data$sr      <- apply(data,1, function(x){sum(x>0)})#species richness
  data$shannon <-diversity(data[,1:length(data)]) # shannon diversity 
  
  data<-data%>%rownames_to_column()
  colnames(data)[1] <- 'sample_id'
  table <- data %>%
    left_join(sample_info,by = 'sample_id')
  return(table)
}



# make colors transparent
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1,], rgb.val[2,], rgb.val[3,],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  return(t.col)
}



# NMDS plot
beta_plot <- function(beta_matrix,sam_table,xlim = c(-0.8,1),ylim = c(-0.4,0.7),yaxt = c(-0.3,0,0.3,0.6),
                    xaxt = c(-1, 0, 1, 2), tag = 'a'){
  par(mar=c(4, 4, 3, 5.9),xpd = T)
  
  
  plot(NA,col= sam_table$col,pch = sam_table$pch,
       xlim = xlim,ylim = ylim,
       ylab = 'NMDS2',xlab = 'NMDS1',
       yaxt="n", xaxt = 'n',font.main = 1,
       cex.axis = 0.7)
  axis(2, at=yaxt,sprintf("%0.1f", yaxt), las=2,tck = -0.03, cex.axis = 0.7)
  axis(1, at=xaxt,sprintf("%0.1f", xaxt), las=1,tck = -0.03, cex.axis = 0.7)
  ordiellipse(beta_matrix,groups=sam_table$status,draw="polygon",col=c('gray50',col_fig),
              label=F,border = NA, alpha= 255*0.55)
  points(beta_matrix$points,col= sam_table$col,pch = 16,cex= 0.8)
  fig_label(tag, font = 2)
  
  legend('topleft', paste('Stress = ', beta_matrix$stress %>% round(3)), bty = 'n', cex = 0.7)
}




# boxplot for diversity and abundance of soil microbes
boxplot_div <- function(table_alpha = alpha_16s, y = 'sr', title = 'Bacteria', ylab = 'Species richness', tag = 'a'){
  table_alpha %>% 
    filter(status != 'back') %>% 
    ggplot() + aes_string(x = 'status', y = y, col = 'status') + geom_boxplot(size = 0.5, outlier.size = 0.5) + 
    scale_color_manual(values = c('black', col_fig)) +
    facet_wrap(~ ster.not,labeller = labeller(ster.not = c('nonster' = 'Live', 'ster' = 'Sterilized'))) + 
    labs(x = NULL, y = ylab, title = title, tag = tag) + 
    theme + theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1, colour = 'black', margin = margin(t = 6)),
                  plot.margin = margin(t = 0.1, r = 0.1, b = 0.3, l = 0.1, unit = "cm"))
}




### figure 4e-g
forest_beta_t <- function(beta_mat = beta_mat_sp_16s_nonster, title = 'Bacteria', 
                          xlim = c(0.38,0.57), xaxt = c(0.4, 0.45, 0.5, 0.55), ylab = F, tag = 'e'){
  col = c(rep('black', 3), rep('gray35',2))
  yaxt = c(3, 2 , 1 , 5.5, 4.5)
  table <- beta_mat %>% 
    group_by(combi_status) %>% 
    summarise(mean = mean(mean_beta), se = sd(mean_beta)/sqrt(n()))
  plot(NA, 
       xlim = xlim,
       ylim = c(0.5, 6), xaxt="n", yaxt = 'n',
       ylab = '', xlab = '', main = title, cex.lab = 0.8)
  segments(table$mean - table$se, yaxt, table$mean + table$se, yaxt, col = col)
  points(table$mean, yaxt, pch =19, col = col, cex = 0.6)
  # y axis
  axis(2, yaxt, labels = rep('', 5),las=1, cex.axis=0.8, las = 2, tck=-0.03)
  if (ylab == T){
    mtext(table$combi_status, side = 2, line = 0.45, at = yaxt, col = col, cex = 0.55, las = 2)
  }
  # x axis
  axis(1, at=xaxt, labels = rep('',length(xaxt)),las=1, cex.axis = 0.8, tck=-0.03)
  mtext(sprintf("%0.2f", xaxt), side = 1, line = 0.6, at = xaxt,  cex = 0.5, las = 1)
  # xlab
  mtext('Soil community dissimilarity', side= 1, line = 2, at=par("usr")[1] + 0.45*diff(par("usr")[1:2]), cex = 0.55)
  # tag
  mtext(tag, side= 3, line= 0.5, at=par("usr")[1] - 0.25*diff(par("usr")[1:2]), cex = 0.7, font = 2)
}



 


# add tag on base plot.
fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}


