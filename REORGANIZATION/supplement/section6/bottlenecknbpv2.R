## Implementation of Algorithm 1
## Main fumction:  minmaxnbp(w), 2 the square weight/distance matrix

###
# functions
###

## Find the path from a vertex v back to s

FindPath_shrinked <- function(v, p, m, s){

	Path = c(v, p[v])
	
	while(1){
		#print(Path)
		Pathtemp = tail(Path, 1)
		if(is.na(Pathtemp)){
			Path = Path[-length(Path)]
		}
		if(Pathtemp==s) break;

		Pathtemp = c(m[Pathtemp], p[m[Pathtemp]])
		Path = c(Path, Pathtemp)
	}
	Path
}
			
## Find M\oplus P, the augmented match			
augmentMatch <- function(.M, .Path){
	while(1){
		e = head(.Path, 2)
		Medge = apply(.M, 1, function(me) prod(me == e) | prod(me== rev(e)))
		if( sum(Medge) > 0){
			.M = .M[!Medge,,drop=FALSE]
		} else {
			.M = rbind(.M, e)
		}
		.Path = .Path[-1]
		if(length(.Path)<=1)
			break		
	}

	.M
}

### Scan vertex u for reassignmet of the uncolored vertices
scan <- function(u, v, p, p_org, l, dplus, dminus, w, V1){
	for(v.prime in names(V1))
		#if(l[v.prime] == 'uncol' & v.prime!=v)
		if(l[v.prime] == 'uncol' & v.prime!=u & v.prime!=v)
		if(max(dplus[u], w[u, v.prime]) <= dminus[v.prime]){
			dminus[v.prime] = max(dplus[u], w[u, v.prime])
			#cat("change: parent of ", v.prime, " to ", u,"\n")
			p[v.prime] = u
			p_org[v.prime] = u
		}

	return(list(p=p, p_org = p_org, dminus = dminus))		
}



## Expand the blossom Bl, starting from v1, to v0
Expand <- function(v1, Bl, v0, B, Pseudo, V_rev, m, p, p_org, s, Bl_s){
	#cat("Expand: v1 = ", v1, ", Bl = ", Bl, ", v0 = ", v0, "\n")
	
	## Special case:
		
	pv1 = p_org[v1]	## This is a node or pseudonode in Bl
	## Find where it is 
	node = V_rev[pv1]
	## Look at the Path in node
	
	B[[node]]	## This can just be B[[Bl]] ## But may not be as well.
	## Find the position of pv1 wrt B[[node]]
	pv1.pos = which(pv1 == B[[node]])
	if(pv1.pos == 1){
		## Need to do something
		Pathtemp = c()
	} else if(pv1.pos %% 2  == 1){ # Go left
		Pathtemp = rev(B[[node]][2:pv1.pos])
	} else {	# Go right
		Pathtemp = B[[node]][pv1.pos:length(B[[node]])]
	} 
	
	Pathtemp = c(Pathtemp, B[[node]][1])
	
	#cat("Pathtemp: ", Pathtemp, "\n")
	
	Pathtemp = ExpandPath(Pathtemp, B, Pseudo, V_rev, m, p, p_org, s, Bl_s)
	## Incomplete
	
	#print(Pathtemp)
	if(tail(Pathtemp, 1) == s){
		return(Pathtemp)
	} else if(is.na(v0)){	## Bl is a blossom containing s
		
		v1 = m[tail(Pathtemp, 1)]
		Pathtemp = c(Pathtemp, m[tail(Pathtemp, 1)])
		#cat("v1: ", v1, "\n")
		Pathtemp1 = FindPath_shrinked(v1, p, m, Bl_s)
		Pathtemp1 = ExpandPath(Pathtemp1, B, Pseudo, V_rev, m, p, p_org, s, Bl_s) 
		Pathtemp = c(Pathtemp, Pathtemp1[-1])
		
	} else if(m[tail(Pathtemp, 1)] == v0){
		return(Pathtemp)
	} else {
		
		v1 = m[tail(Pathtemp, 1)]
		Pathtemp = c(Pathtemp, m[tail(Pathtemp, 1)])
		#cat("Now, v1: ", v1, "\n")
		Pathtemp1 = FindPath_shrinked(v1, p, m, p[v0])
		Pathtemp1 = ExpandPath(Pathtemp1, B, Pseudo, V_rev, m, p, p_org, s, Bl_s) 
		Pathtemp1 = Pathtemp1[1:(which(Pathtemp1==v0)-1)]
		Pathtemp = c(Pathtemp, Pathtemp1[-1])
	}
	return(Pathtemp)
	
}

## Expand a Path, Expand and ExpandPath interacts with each other.
ExpandPath <- function(Path, B, Pseudo, V_rev, m, p, p_org, s, Bl_s){

	#print(Path)
	L = length(Path)
	
	if(L<1) return(Path)
	
	Path_new = c()
	for(i in 1:length(Path)){
		if(Path[i] == s){
			Path_new = c(Path_new, Path[i])
			break;
		}
		if(Pseudo[Path[i]] != 1)
			Path_new = c(Path_new, Path[i])
		if(Pseudo[Path[i]] == 1){
			if(!is.na(m[Path[i]]) & m[Path[i]] == Path[i-1]){
				Pathtemp = B[[Path[i]]][1]
				while(1){
					if(Pseudo[Pathtemp] == 1){
						Pathtemp = B[[Pathtemp]][1]
					} else break;
				}
			} else if(i==length(Path)){
				#Path_new = c(Path_new, Path[i])
				Pathtemp = Expand(Path[i-1], Path[i], m[Path[i]], B, Pseudo, V_rev, m, p, p_org, s, Bl_s)
			} else {
				Pathtemp = Expand(Path[i-1], Path[i], Path[i+1], B, Pseudo, V_rev, m, p, p_org, s, Bl_s)
			}
			Path_new = c(Path_new, Pathtemp)
		}
	}
	
	return(Path_new)
}
	
	

	
	
#######################


minmaxnbp <- function(w){

	n = dim(w)[1]
	V = as.character(1:n)
	names(V) = V
	rownames(w) = colnames(w) = V
	
	temp = Sys.time()

	M = matrix(NA, 0, 2)

	count = 0

	while(1){

		## search s
		s = 1
		while(s<=n){
			if(!(V[s] %in% as.vector(M))){
				break;
			}
			s = s+1
		}

		s = V[s] #as.character(s)

		s

		#s = "99"
		###############################################
		## Initialize
		dplus = rep(NA, n)
		dminus = rep(NA, n)
		names(dplus) = names(dminus) = V

		## l: store color: 'col', 'uncol' or 'inline'
		## p: store potential parent
		## m: matched vertex
		l <- rep("uncol", n)
		p <- rep(NA, n)	
		names(l) = names(p) = V
		m <- rep(NA, n)
		names(m) = V
		
		w1 = w
		rownames(w1) = colnames(w1) = V
		V1 = as.character(1:n)
		names(V1) = V
		V_rev = V
		V = V1
		Pseudo = rep(0, n)
		names(Pseudo) = V
		
		
		B = list()
		Bmaxvalues = list()
		MA_matches = M
		
		l[s] = "col"
		dminus[which(l!='col')] = w[s, which(l!='col')]
		dplus[s] = Inf
		dminus[s] = -Inf
		p[which(l!='col')] = s	## everyone is connected to s
		## p_org is the original parent of the vertex. 
		p_org = p
		
		Bl_s = s
		
		###############################################
		## Pick the best

		while(1){
			#print(dminus)
			#print(dplus)
			delta1 = min(dminus[l[names(V1)]=='uncol'], na.rm=TRUE)
			delta2 = min(pmax(dplus, dminus)[l[names(V1)]=='col'], na.rm=TRUE)
			#delta2 = min((dplus)[l[names(V1)]=='col'], na.rm=TRUE)

			## Will add two 
			if(delta1 < delta2)
				v = Find(function(x){ if(l[x]!='uncol') return(FALSE)
									x = dminus[x]; ifelse(is.na(x), FALSE, x==delta1)},
								names(V1))
				#which(dminus == delta1)
			if(delta2 <= delta1)
				v = Find(function(x){ if(l[x]!='col') return(FALSE)
								x = max(dplus[x], dminus[x]); 
								ifelse(is.na(x), FALSE, x==delta2)}, names(V1))
				#v = Find(function(x){ x = dplus[x]; ifelse(is.na(x), FALSE, x==delta2)}, names(V1))
				#which.min(dplus)

			v
			delta2 <= delta1
			#cat("delta1: ", delta1, "delta2: ", delta2,"\n")
			## If delta1 is smaller?

			if(delta1 < delta2){
				#MA_matches = MA(M, B)

				if(!(v %in% MA_matches)){
					#print("Found one!")
					
					#print(M[order(as.numeric(M[,1])),])	
					## Find the shortest augmenting path by 
					## expanding the path from v to s
					
					shortestPath = FindPath_shrinked(v, p, m, s=Bl_s)
					shortestPath = ExpandPath(shortestPath, B, Pseudo, V_rev, m, p, p_org, s, Bl_s)
					
					#print(shortestPath)
					## Augment the match with the path to create the augmented match
					.M = augmentMatch(M, shortestPath)
					
					#print(.M[order(as.numeric(.M[,1])),])
					break;					

				} else {
					
					u = union(MA_matches[MA_matches[,1]==v,2,drop=FALSE], 
								MA_matches[MA_matches[,2]==v,1,drop=FALSE])

					#u = as.character(u)
					m[u] = v
					l[u] = "col"
					l[v] = 'inline'
					#p_org[u] = v
					#dminus[u] = -Inf

					#dminus[v] = Inf
					dplus[u] = delta1
										
					#cat(p[v], v, u,"\n")
					#cat(p_org[v], v, u,"\n")
					
					## Scan u, for reassignment of uncolored vertices 
					res = scan(u, v, p, p_org, l, dplus, dminus, w1, V1)					
					dminus = res$dminus
					p = res$p
					p_org = res$p_org

					#delta2=delta1
					#cat("Visit: ")
					#cat(p[v], v, u,"\n")
					
					
				}
			} else { #cat('--------------------------------\n\n')
					#print("Create a blossom!") 
					#break;
					#if(length(B)>1)
					#	break;
					
					# Find blossom between v and p[v], using m and p
					#print(c(p[v], v))
					
					Path1 = c(v, FindPath_shrinked(m[v], p, m, s=Bl_s))	#V_rev[s]))
					if(p[v] == V_rev[s] | is.na(m[p[v]])){
						Path2 = c(p[v])#, FindPath_shrinked(m[p[v]], p, m, s=V_rev[s]))
					} else {
						Path2 = c(p[v], FindPath_shrinked(m[p[v]], p, m, s=Bl_s))
					}
					
					topnode.pos1 = which(Path1 %in% Path2)
					topnode.pos2 = which(Path2 %in% Path1) 
					
					b = Path1[1:topnode.pos1[1]]
					if(topnode.pos2[1]>1)
						b = c(rev(Path2[1:(topnode.pos2[1]-1)]), b)
					#if(topnode.pos2[1]==1)					
					b = rev(b)	
					## b is the blossom		
					#print(b)
					
					## 	Attach blossom to list B,
					#cat("\n")
					newNodename = paste0('B',length(B)+1)
					B[[newNodename]] = b
					
					
					## b is now a pseudo-node.
					## Define the parent and the matched vertex of the 
					##	pseudo-node.
					#m[newNodename] = m[p[v]]
					m[newNodename] = m[b[1]]
					
					if(is.na(m[newNodename]))
						Bl_s = newNodename
						
					#p[newNodename] = p[p[v]]
					#p[newNodename] = p[m[p[v]]]
					p[newNodename] = p[m[b[1]]]
					p[p %in% b] = newNodename
					
					## Shrink the pseudo-node, b
					keepNodes = setdiff(names(V1), b)
					V1 = V1[keepNodes]
					V1[newNodename] = paste(b, collapse=',')
					
					Pseudo[newNodename] = 1
					
					## Store the original information about the 
					## psedo-node for reverse search
					V_rev[b] = newNodename

					## Color the pseudo-node as 'colored'
					l[newNodename] = 'col'
					
			
					## Assign the cost of the pseudo nodes.
					dplus = dplus[keepNodes]
					dplus[newNodename] = delta2
				
					dminus = dminus[keepNodes]
					dminus[newNodename] = Inf
					
					## Scan the pseudo nodes
					#cat("dplus[newNodename]: ", dplus[newNodename],"\n")
					#print(names(V1))
					for(v.prime in names(V1)){
						if(l[v.prime] == 'uncol'){
							dtemp = sapply(intersect(V, b),
									function(u) 
										max(dplus[newNodename], w[u, v.prime]))
							names(dtemp) = intersect(V, b)
							
							if(min(dtemp) < dminus[v.prime]){
								dminus[v.prime] = min(dtemp)
								#cat("changed: ", v.prime, ": ", p_org[v.prime])
								p_org[v.prime] = names(which.min(dtemp))
								#cat(", ", p_org[v.prime], "\n")
								p[v.prime] = newNodename
							}
						}
					}
					
					#cat("\n\n")
					
					#if(length(B)>0)
					#	break;
					
					}
		}	
		

		if(sum(table(.M)>1)) {
			print("Error: Unit multiple times in the match!")
			break;
		}
		#if(count == 49) break;
		
		M = .M
		count = count + 1
		dist_max = apply(M, 1, function(b) w[b[1], b[2]])

	
		#cat('nrow(M): ', nrow(M),'\n')
		#cat(count, ": ", max(dist_max), ": Next\n\n")


		#if(count == 15) break;

		if(nrow(M) == n/2)
			break;

	}

	Sys.time() - temp

	M
}
