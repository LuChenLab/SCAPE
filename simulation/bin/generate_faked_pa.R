library(dplyr)
args <- commandArgs(T)
input_bed <- args[1]
output_bed <- args[2]
gap_mean <- args[3]

df <- read.table(input_bed, sep = '\t', stringsAsFactor=F)

df %>% group_by(V7) %>% top_n(1, V5) -> df




df_lst <- split(df, df$V7)
# sample_ind <- sample(1:length(df_lst),100)
# df_lst <- df_lst[sample_ind]


gap_mean <- as.numeric(gap_mean)
gap_sd <- 10

df_lst <- lapply(df_lst, function(x){
	offset = 0
	tmp <- x
	times <- sample(1:4,1)
	for (time in 1:times) { 
		new_line = x

		offset <- offset + gap_mean + round(gap_sd * runif(1))
		tmp_name <- unlist(strsplit(x$V4, split = '_'))
		tmp_name[1] <- paste(tmp_name[1], glue::glue('fake{time}'), sep = '')
		tmp_name <- paste(unlist(tmp_name), collapse='_')
		new_line$V4 <- tmp_name

		if (x$V6 == "+") {

			new_line$V3 <- new_line$V3 - offset
			if (new_line$V3 <= x$V2) {
				break
			}

			} else {

			new_line$V2 <- new_line$V2 + offset
			if (new_line$V2 >= x$V3) {
				break
			}

			}
		tmp <- rbind(tmp, new_line)
		}

	return(tmp)
	})

df <- do.call(rbind, df_lst)
write.table(df, file = output_bed, sep = '\t', quote=F, row.names=F, col.names=F)

