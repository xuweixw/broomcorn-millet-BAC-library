parsePool <- function(x_axis, y_axis, dim = 48){
  # the arragement of 6 384-plates,generating 48 X 2 pools
  #   1 2
  #   3 4
  #   5 6
  # writen in July 25,2017
  well_L <- c()
  a <- LETTERS[1:16]
  for (k in c(1,3,5)){
    for(i in a){
      for (j in c("01","02","03","04","05","06","07","08","09",10:24)){
        well_L <- c(well_L, paste0(k,i,j))
      }
    }
  }
  plate_L <- matrix(well_L,ncol=24,byrow = TRUE)
  
  well_R <-c()
  for (z in c(2,4,6)){
    for(x in a){
      for (y in c("01","02","03","04","05","06","07","08","09",10:24)){
        well_R <- c(well_R, paste0(z,x,y))
      }
    }
  }
  plate_R <- matrix(well_R,ncol=24,byrow = TRUE)
  
  plate_48P <- data.frame(plate_L, plate_R)
  colnames(plate_48P) <- paste0(1:48)
  return(as.vector(plate_48P[x_axis,y_axis]))
}


X <- c("03","11","19","27","35","43","04","13","21","29","36","44")
Y <- c("03","11","19","27","35","43","04","12","20","28","36","44")


barcode <- function(num){
	# num means the number of 384-plate

	# One 384-plate
	row <- LETTERS[1:16]
	col <- c("01","02","03","04","05","06","07","08","09",10:24)
	id_element <- c()
	for (row_letter in row){
		for(col_letter in col){
			id_element <- c(id_element,paste(row_letter, col_letter, sep=""))
        }
    }

    # lots of 384-plate
    id <- c()
    for (n in 1:num){
    	for(id_e in id_element){
    		id <- c(id, paste(n, id_e, sep=""))
    	}
    }

    return(id)
}
