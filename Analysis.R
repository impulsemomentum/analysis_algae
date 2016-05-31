
library(car)
library(DAAG)
library(DMwR)
library(lattice)
draw <- function(list,list1,file,title)
{
	ppath = file
	png(file=ppath,width=1200,height=1000)
	split.screen(c(3,2))
	par(pin=c(4,4))
	screen(1) 
	hist(list,main = paste(title," ",sep =""))
	screen(3)
	qqPlot(list,main=paste(title," ",sep =""))
	screen(5)
	boxplot(list, main=paste(title," ",sep =""), na.rm = TRUE)
	
	screen(2) 
	hist(list1,main = paste(title," original",sep =""))
	screen(4)
	qqPlot(list1,main=paste(title," original",sep =""))
	screen(6)
	boxplot(list1, main=paste(title," original",sep =""), na.rm = TRUE)
	dev.off()
	
}

draw_all <- function(Adata,data,dir)
{
	draw(Adata$mxPH,data$mxPH,paste(dir,"mxPH.png",sep =""),"mxPH")
	draw(Adata$mnO2,data$mnO2,paste(dir,"mnO2.png",sep =""),"mnO2")
	draw(Adata$Cl,data$Cl,paste(dir,"Cl.png",sep =""),"Cl")
	draw(Adata$NO3,data$NO3,paste(dir,"NO3.png",sep =""),"NO3")
	draw(Adata$NH4,data$NH4,paste(dir,"NH4.png",sep =""),"NH4")
	draw(Adata$oPO4,data$oPO4,paste(dir,"oPO4.png",sep =""),"oPO4")
	draw(Adata$PO4,data$PO4,paste(dir,"PO4.png",sep =""),"PO4")
	draw(Adata$Chla,data$Chla,paste(dir,"Chla.png",sep =""),"Chla")

	draw(Adata$a1,data$a1,paste(dir,"a1.png",sep =""),"a1")
	draw(Adata$a2,data$a2,paste(dir,"a2.png",sep =""),"a2")
	draw(Adata$a3,data$a3,paste(dir,"a3.png",sep =""),"a3")
	draw(Adata$a4,data$a4,paste(dir,"a4.png",sep =""),"a4")
	draw(Adata$a5,data$a5,paste(dir,"a5.png",sep =""),"a5")
	draw(Adata$a6,data$a6,paste(dir,"a6.png",sep =""),"a6")
	draw(Adata$a7,data$a7,paste(dir,"a7.png",sep =""),"a7")
}


clean_NA <-function(Adata)
{
	Adata[!complete.cases(Adata),]
	Adata <- na.omit(Adata)
	Adata
}

central_NA <- function(Adata)
{
	#data(Adata)
	Adata <- Adata[-manyNAs(Adata)]
	Adata <- centralImputation(Adata)
	Adata
}

knn_NA <- function(Adata)
{
	Adata <- knnImputation(Adata,10)
}

relation_NA <- function(Adata)
{
	symnum(cor(Adata[,4:18],use='complete.obs'))
	lm(formula=PO4~oPO4, data=Adata)
	Adata = Adata[-manyNAs(Adata),]
	Adata
}


location = "Analysis.txt"

data <- read.table(location,col.name = c('season','size','speed','mxPH','mnO2','Cl','NO3','NH4','oPO4','PO4','Chla','a1','a2','a3','a4','a5','a6','a7'),na.string=('XXXXXXX'))

summary(data)


clean_NA_data = clean_NA(data)
dir.create("clean_NA")
draw_all(clean_NA_data,data,"clean_NA/")
write.csv(clean_NA_data, file = "clean_NA/OmitedData.csv")

dirs = "clean_NA/"
png(file=paste(dirs,"a_1.png",sep=""),width=600,height=500)
bwplot(a1~size,data=data,ylab='River Size',xlab='a1')
png(file=paste(dirs,"a_7.png",sep=""),width=600,height=500)
bwplot(a7~size,data=data,ylab='River Size',xlab='a7')
png(file=paste(dirs,"a_2.png",sep=""),width=600,height=500)
bwplot(a2~size,data=data,ylab='River Size',xlab='a2')
png(file=paste(dirs,"a_3.png",sep=""),width=600,height=500)
bwplot(a3~size,data=data,ylab='River Size',xlab='a3')
png(file=paste(dirs,"a_4.png",sep=""),width=600,height=500)
bwplot(a4~size,data=data,ylab='River Size',xlab='a4')
png(file=paste(dirs,"a_5.png",sep=""),width=600,height=500)
bwplot(a5~size,data=data,ylab='River Size',xlab='a5')
png(file=paste(dirs,"a_6.png",sep=""),width=600,height=500)
bwplot(a6~size,data=data,ylab='River Size',xlab='a6')


central_NA_data = central_NA(data)
dir.create("central_NA")
draw_all(central_NA_data,data,"central_NA/")
write.csv(central_NA_data, file = "central_NA/CentralImputationData.csv")
dirs = "central_NA/"
png(file=paste(dirs,"a_1.png",sep=""),width=600,height=500)
bwplot(a1~size,data=data,ylab='River Size',xlab='a1')
png(file=paste(dirs,"a_7.png",sep=""),width=600,height=500)
bwplot(a7~size,data=data,ylab='River Size',xlab='a7')
png(file=paste(dirs,"a_2.png",sep=""),width=600,height=500)
bwplot(a2~size,data=data,ylab='River Size',xlab='a2')
png(file=paste(dirs,"a_3.png",sep=""),width=600,height=500)
bwplot(a3~size,data=data,ylab='River Size',xlab='a3')
png(file=paste(dirs,"a_4.png",sep=""),width=600,height=500)
bwplot(a4~size,data=data,ylab='River Size',xlab='a4')
png(file=paste(dirs,"a_5.png",sep=""),width=600,height=500)
bwplot(a5~size,data=data,ylab='River Size',xlab='a5')
png(file=paste(dirs,"a_6.png",sep=""),width=600,height=500)
bwplot(a6~size,data=data,ylab='River Size',xlab='a6')


knn_NA_data = knn_NA(data)
dir.create("knn_NA")
draw_all(knn_NA_data,data,"knn_NA/")
write.csv(knn_NA_data, file = "knn_NA/knnImputationData.csv")
dirs = "knn_NA/"
png(file=paste(dirs,"a_1.png",sep=""),width=600,height=500)
bwplot(a1~size,data=data,ylab='River Size',xlab='a1')
png(file=paste(dirs,"a_7.png",sep=""),width=600,height=500)
bwplot(a7~size,data=data,ylab='River Size',xlab='a7')
png(file=paste(dirs,"a_2.png",sep=""),width=600,height=500)
bwplot(a2~size,data=data,ylab='River Size',xlab='a2')
png(file=paste(dirs,"a_3.png",sep=""),width=600,height=500)
bwplot(a3~size,data=data,ylab='River Size',xlab='a3')
png(file=paste(dirs,"a_4.png",sep=""),width=600,height=500)
bwplot(a4~size,data=data,ylab='River Size',xlab='a4')
png(file=paste(dirs,"a_5.png",sep=""),width=600,height=500)
bwplot(a5~size,data=data,ylab='River Size',xlab='a5')
png(file=paste(dirs,"a_6.png",sep=""),width=600,height=500)
bwplot(a6~size,data=data,ylab='River Size',xlab='a6')
dev.off()


relation_NA = relation_NA(data)
dir.create("relation_NA")
draw_all(relation_NA,data,"relation_NA/")
write.csv(relation_NA, file = "relation_NA/linearDefaultData.csv")
dirs = "relation_NA/"
png(file=paste(dirs,"a_1.png",sep=""),width=600,height=500)
bwplot(a1~size,data=data,ylab='River Size',xlab='a1')
png(file=paste(dirs,"a_7.png",sep=""),width=600,height=500)
bwplot(a7~size,data=data,ylab='River Size',xlab='a7')
png(file=paste(dirs,"a_2.png",sep=""),width=600,height=500)
bwplot(a2~size,data=data,ylab='River Size',xlab='a2')
png(file=paste(dirs,"a_3.png",sep=""),width=600,height=500)
bwplot(a3~size,data=data,ylab='River Size',xlab='a3')
png(file=paste(dirs,"a_4.png",sep=""),width=600,height=500)
bwplot(a4~size,data=data,ylab='River Size',xlab='a4')
png(file=paste(dirs,"a_5.png",sep=""),width=600,height=500)
bwplot(a5~size,data=data,ylab='River Size',xlab='a5')
png(file=paste(dirs,"a_6.png",sep=""),width=600,height=500)
bwplot(a6~size,data=data,ylab='River Size',xlab='a6')
dev.off()
