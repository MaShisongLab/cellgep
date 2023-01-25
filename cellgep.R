cellgep <- function ( expression, neighbor, selected_gep, dynamic.cutoff=T)
{
	require(graph)
	require(RBGL)
	
	idx_gep = 1:dim(expression)[1] %in% selected_gep
	b = expression
	b[!idx_gep,] = 0
	b_max = apply(b,1,max)
	b_mean = apply(b,1,mean)
	b_sd = apply(b,1,sd)
	c = b
	geps_to_test = selected_gep
	
	if ( dynamic.cutoff){
		cutoff_info_colnames = c('GEP','Level','Cutoff','Found','CellNum','Num1stCC','Num2ndCC','Fraction1stCC','Fraction1&2CC','NumCC')
		cat('GEP\tLevel\tCutoff\tFound\tCellNum\tNum1stCC\tNum2ndCC\tFraction1stCC\tFraction1&2CC\tNumCC\n')
	}else{
		cutoff_info_colnames = c('GEP','Cutoff')
	}

	for( i in geps_to_test){
		if ( dynamic.cutoff ){
			cutoff_result = get_dynamic_cutoff(b, i, neighbor, b_max, b_mean, b_sd)
			cutoffdata = cutoff_result$result
			cutoff = as.numeric(cutoffdata[1])
			cutoffdata[1] = round(cutoff,3)
			cat(row.names(expression)[i], round(cutoff / (b_max[i] + 0.0000000001), 2), cutoffdata,'\n',sep='\t')
			cutoffdata = c(row.names(expression)[i],round(cutoff / (b_max[i] + 0.0000000001), 2),cutoffdata)
			
			int.results = cutoff_result$int.results
			int.results = cbind(rep(row.names(expression)[i],dim(int.results)[1]),int.results)
			
			if ( cutoff_result$multiple.rows == "N"){
				int.results = int.results[1,]
			}

			if( exists('intermediate.results')){
				intermediate.results = rbind(intermediate.results, int.results)
			}else{
				intermediate.results = int.results
			}
		}else{
			cutoff = max(0.3, b_max[i]/2)
			cutoffdata = c(row.names(expression)[i],cutoff)
		}

		if( exists('cutoff_info')){
			cutoff_info = rbind(cutoff_info, cutoffdata)
		}else{
			cutoff_info = cutoffdata
		}

		idx = b[i,] < cutoff
		c[i, idx] = 0
	}
	max_identity = apply(c,2,which.max)
	max_value = apply(c,2,max)
	max_identity[max_value == 0] = "na"

	marker_gep = max_identity

	nn = neighbor
	neighbor_row_num = dim(neighbor)[1]
	neighbor_col_num = dim(neighbor)[2]
	for ( i in 1:neighbor_row_num){
		nn[i,] = marker_gep[neighbor[i,]]
	}

	new_id <- as.character(nn[,1])
	is_original_label <- new_id
	cell_expression <- 1:neighbor_row_num
	cell_gep_mean <- 1:neighbor_row_num
	cell_gep_sd <- 1:neighbor_row_num

	for ( i in 1:neighbor_row_num){
		a = as.character(nn[i,1])
		b = as.character(nn[i,2:neighbor_col_num])
		if ( a %in% b & !(a == "na") ){     # & table(as.character(b))[a] > 0
			is_original_label[i] = "Y"
			cell_expression[i] = round(expression[as.numeric(a),i],3)
			cell_gep_mean[i] <- round(b_mean[as.numeric(a)],3)
			cell_gep_sd[i] <- round(b_sd[as.numeric(a)],3)
		}else{
			b1 = 1:length(b)
			b2 = cbind(b1,b1)
			b2[,1] = 1
			t = aggregate(b2,by=list(b),sum)
			idx = t[,1] != "na"
			t = t[idx,]
			t = t[order(-t[,2],t[,3]),]
			if( sum(idx) > 0 & t[1,2] > 2){
				id_m = as.character(t[1,1])
				cell_exp = expression[as.numeric(id_m),i]
				if (cell_exp >= 0.3 & cell_exp >= b_max[as.numeric(id_m)]/100 ) #b_mean[as.numeric(id_m)] + b_sd[as.numeric(id_m)])
				{
					new_id[i] = id_m
					cell_expression[i] = round(expression[as.numeric(id_m),i],3)
					cell_gep_mean[i] <- round(b_mean[as.numeric(id_m)],3)
					cell_gep_sd[i] <- round(b_sd[as.numeric(id_m)],3)
				}else{
					new_id[i] = "na"
					cell_expression[i] = "na"
					cell_gep_mean[i] <- "na"
					cell_gep_sd[i] <- "na"
				}
			}else{
				new_id[i] = "na"
				cell_expression[i] = "na"
				cell_gep_mean[i] <- "na"
				cell_gep_sd[i] <- "na"
			}
			if ( new_id[i] != a) 
			{
				is_original_label[i] = "N"
			}else{
				is_original_label[i] = "Y"
			}
		}
	}

	# Two rounds of annotating 'na' cells by neighbor-voting
	for ( round in 1:2){
		for ( i in 1:neighbor_row_num){
			if (new_id[i] == "na"){
				b = as.character(new_id[as.numeric(neighbor[i,2:neighbor_col_num])])
				b1 = 1:length(b)
				b2 = cbind(b1,b1)
				b2[,1] = 1
				t = aggregate(b2,by=list(b),sum)
				idx = t[,1] != "na"
				t = t[idx,]
				t = t[order(-t[,2],t[,3]),]
				if( sum(idx) > 0 & t[1,2] > 2){
					id_m = as.character(t[1,1])
					cell_exp = expression[as.numeric(id_m),i]
					if (cell_exp >= 0.3 & cell_exp >= b_max[as.numeric(id_m)]/100)
					{
						new_id[i] = id_m
						cell_expression[i] = round(expression[as.numeric(id_m),i],3)
						cell_gep_mean[i] <- round(b_mean[as.numeric(id_m)],3)
						cell_gep_sd[i] <- round(b_sd[as.numeric(id_m)],3)
						is_original_label[i] = "N"
					}
				}
			}	
		}
	}

	idx = new_id != "na"
	gep = new_id
	gep[idx] = row.names(expression)[as.numeric(new_id[idx])]

	geps_expressed = new_id
	for ( i in 1:neighbor_row_num){
		idx = which(c[,i] > 0)
		if( length(idx) > 0){
			aa = expression[idx,i]
			names(aa) = paste(row.names(expression)[idx],round(expression[idx,i],2),sep="-")
			aa = sort(aa,decreasing=T)
			aa = paste(names(aa),collapse=", ")
		}else{
			aa = "na"
		}
		geps_expressed[i] = aa
	}

	e = cbind(gep,new_id,is_original_label,cell_expression,cell_gep_mean,cell_gep_sd,geps_expressed,nn)
	colnames(e)[1:8] = c("GEP","GEP_id","is_original","Expression","GEP_mean","GEP_sd","GEPs_expressed","GEP_Ori")
	colnames(e)[9:(neighbor_col_num + 7)] = paste("N",1:(neighbor_col_num - 1),sep="")
	e = as.data.frame(e)

	colnames(cutoff_info) = cutoff_info_colnames
	row.names(cutoff_info) = NULL
	cutoff_info = as.data.frame(cutoff_info)
	if (dynamic.cutoff){
		if( length(intermediate.results > 9)){
			colnames(intermediate.results) = c('GEP','Cutoff','Qulified','NumCell','Num1stCC','Num2ndCC','Fraction1stCC','Fraction1&2CC','NumCC')
			row.names(intermediate.results) = NULL
		}else{
			names(intermediate.results) = c('GEP','Cutoff','Qulified','NumCell','Num1stCC','Num2ndCC','Fraction1stCC','Fraction1&2CC','NumCC')
		}
		intermediate.results = as.data.frame(intermediate.results)
		intermediate.results[,2] = round(as.numeric(intermediate.results[,2]),3)
		intermediate.results[,4] = as.numeric(intermediate.results[,4])
		intermediate.results[,5] = as.numeric(intermediate.results[,5])
		intermediate.results[,6] = as.numeric(intermediate.results[,6])
		intermediate.results[,7] = as.numeric(intermediate.results[,7])
		intermediate.results[,8] = as.numeric(intermediate.results[,8])
		intermediate.results[,9] = as.numeric(intermediate.results[,9])
	}else{
		intermediate.results = NULL
	}
		
	return(list(annotation=e,cutoff=cutoff_info,intermediate.results = intermediate.results,expression.filtered = c))
}

get_dynamic_cutoff <- function(b, i, neighbor, b_max, b_mean, b_sd)
{
	require(graph)
	require(RBGL)
	
	c = c(seq(0.95,0.35, by= -.05) * b_max[i],  b_max[i] / 3)
	c = sort(unique(c), decreasing = T)
	low = max(0.3, b_max[i] / 2)
	c = c[c > low]
	c = c(c,low)

	skip = 0
	for ( k in 1:length(c)){
		if (skip == 0){
			ccc = k
			cutoff = c[k]
			idx = b[i,] >= cutoff
			n = neighbor[idx,]
			if( sum(idx) > 1){
				N <- as.character(n[,1])
				gR <- new("graphNEL", nodes = N)
				edges = n[,1:2]
				for ( j in 3:dim(n)[2]){
					edges <- rbind(edges,n[,c(1, j)])
				}
				i2 = as.character(edges[,1]) %in% N & as.character(edges[,2]) %in% N
				if( sum(i2) > 1 ){
					edges = edges[ i2,]
					gX <- addEdge(as.character(edges[,1]),as.character(edges[,2]),gR)
					cc <- connectedComp(gX)
					size = 0
					for ( j in 1:length(cc)) { size = c(size, length(cc[[j]])) }
					N = sum(idx)
					size = sort(size, decreasing=T)

					if ( (size[1] / N >= 0.98 & size[1] >= 5) | (sum(size[1:2]) == N & size[2] >= 10 )){
						qualified = "Y"
					}else{
						qualified = "N"
					}

					tt = c(cutoff, qualified, N, size[1], size[2], round( size[1] / N, 3), round( sum(size[1:2])/N, 3), length(cc))
				}else{
					tt = c(cutoff, "N", sum(idx), 0, 0, 0, 0, 0)
				}
			}else{
				tt = c(cutoff, "N", sum(idx), 0, 0, 0, 0, 0)
			}
			
			if ( k == 1 ){
				ttt = tt
			}else{
				ttt = rbind(ttt,tt)
			}

			if ( (cutoff <= b_max[i] * 0.5) & sum(idx) >= 1000){
				skip = 1;
			}
		}
	}

	if( length(c) > 1){
		if ( "Y" %in% ttt[,2]){
			if( ttt[ccc,2] == "N"){
				skip = 0;
				for( k in (ccc - 1):1){
					if ( skip == 0){
						if ( ttt[k,2] == "Y"){
							result = ttt[k,]
							skip = 1
						}
					}
				}
			}else{
				skip = 0;
				if( ccc < length(c)){
					for( k in (ccc + 1):length(c)){
						if ( skip == 0){
							cutoff = c[k]
							idx = b[i,] >= cutoff
							n = neighbor[idx,]
							N <- as.character(n[,1])
							gR <- new("graphNEL", nodes = N)
							edges = n[,1:2]
							for ( j in 3:dim(n)[2]){
								edges <- rbind(edges,n[,c(1, j)])
							}
							i2 = as.character(edges[,1]) %in% N & as.character(edges[,2]) %in% N
				
							edges = edges[ i2,]
							gX <- addEdge(as.character(edges[,1]),as.character(edges[,2]),gR)
							cc <- connectedComp(gX)
							size = 0
							for ( j in 1:length(cc)) { size = c(size, length(cc[[j]])) }
							N = sum(idx)
							size = sort(size, decreasing=T)

							if ( (size[1] / N >= 0.98 & size[1] >= 5) | (sum(size[1:2]) == N & size[2] >= 10)){
								qualified = "Y"
							}else{
								qualified = "N"
								skip = 1
							}

							tt = c(cutoff, qualified, N, size[1], size[2], round( size[1] / N, 3), round( sum(size[1:2])/N, 3), length(cc))
							ttt = rbind(ttt,tt)

							if ( cutoff <= 0.3 | (cutoff <= (b_mean[i] + b_sd[i]) & (N >= 10000 & N / dim(neighbor)[1] >= 0.1))) {
								#skip = 1
							}
						}
					}
				}
	
				skip = 0
				for ( k in (dim(ttt)[1]):1){
					if (skip == 0){
						if ( ttt[k,2] == "Y")
						{
							result = ttt[k,]
							skip = 1
						}
					}
				}
			}
		}else{
		 	result = c(b_max[i] + 0.001, "N", 0, 0, 0, 0, 0, 0)
		}
		multiple.rows = "Y"
	}else{
		result = c(b_max[i] + 0.001, "N", 0, 0, 0, 0, 0, 0)
		tt = c(b_max[i], "N", 0, 0, 0, 0, 0, 0)
		ttt = rbind(result, tt)
		multiple.rows = "N"
	}

	return(list(result = result, multiple.rows = multiple.rows, int.results = ttt))
}

average_by_gep <-
function( expression, gep, gep_origin = "na", expression_origin = "na")
{
	gene = row.names(expression)
	
	gep_input = gep
	idx = gep[,1] %in% gene
	gep = gep[idx,]
	
	t = table(gep[,2])
	gene_count = t

	t = t[t>=10]
	gep_filtered = gep
	
	r = gene %in% gep[gep[,2] == names(t)[1],1]
	n = sum(r)
	if ( n > 0) { r = r / n }
	a = r
	if( length(names(t)) > 1){
		for( i in 2:length(names(t))){
			r = gene %in% gep[gep[,2] == names(t)[i],1]
			n = sum(r)
			if ( n > 0 ) { r = r / n }
			a = rbind( a, r)
		}
		colnames(a) = gene
		row.names(a) = names(t)
	}else{
		a = rbind(a,a)
		colnames(a) = gene
		row.names(a)[1] = names(t)[1]
		row.names(a)[2] = names(t)[1]
		a = a[1,]
	}
	b = apply(a,1,function(x) sum(x>0))
	if ( sum(b>=5) >= 1){
		a = a[b>=5,]
		require(edgeR)
		for( i in 0:floor((dim(expression)[2] - 1) / 10000)){
			aa = i*10000 + 1
			bb = aa + 9999
			bb = min(bb,dim(expression)[2])
			cc = expression[,aa:bb]
			cc = cpm(cc)/100
			cc = log(cc + 1, 2)
			e = a %*% cc
			if ( i == 0){
				all = e
			}else{
				all = cbind(all,e)
			}
			cat(i, '\n')
		}

	}else{
		cat('No mdoule found.\n')
		all = 0
	}
	colnames(all) = colnames(expression)

	return(list(expression = all, gep_input = gep_input, gep_filtered = gep_filtered, gene_count = gene_count, gep_origin = gep_origin, expression_origin = expression_origin))
}

