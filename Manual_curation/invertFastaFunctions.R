getContigs = function(chr, general, contigs){

        coordinates = general[general[,1] == chr & general[,4] %in% contigs,]

        return(coordinates)

}

getDoubleContigs = function(chr, general, double_contigs){

        # deal with contigs in tandem that must be flipped together
        double_contigs = strsplit(x = double_contigs, split = ",")
        double_coordinates = data.frame()

        # begin for loop
        for(i in 1:length(double_contigs)){

                begin = general[general[,1] == chr & general[,4] %in% double_contigs[[i]][1],2]
                end = general[general[,1] == chr & general[,4] %in% double_contigs[[i]][length(double_contigs)],3]
                double_coordinates = rbind(double_coordinates, data.frame(V1 = chr, V2 = begin, V3 = end, V4 = paste(double_contigs[[i]], collapse = ",")))

        }
        #end for loop

        return(double_coordinates)

}

getDoubleCoordinates = function(chr, general, double_contigs, multiple){

        double_contigs = unlist(strsplit(x = multiple, split = ","))
        begin = general[general[,1] == chr & general[,4] %in% double_contigs[1],2]
        end = general[general[,1] == chr & general[,4] %in% double_contigs[length(double_contigs)],3]

        double_coordinates = data.frame(V1 = chr, V2 = begin, V3 = end, V4 = paste(double_contigs, collapse = ","))

        return(double_coordinates)
}

invertCoordinates = function(double_coordinates){

        new_coordinates = double_coordinates
        min = min(new_coordinates[,2])
        cumulative_length = 0

        # begin for loop
        for(i in 1:nrow(new_coordinates)){

                length_itself = abs(double_coordinates[i,3] - double_coordinates[i,2])
                cumulative_length = abs(double_coordinates[i,3] - double_coordinates[nrow(double_coordinates), 3])
                new_coordinates[i,2] = min + length_itself + cumulative_length
                new_coordinates[i,3] = min + cumulative_length

        }
        # end for loop

        new_coordinates = new_coordinates[(nrow(new_coordinates):1),c(1,3,2,4)]

        return(new_coordinates)

}

getMultipleContigs = function(chr, general, multiple_contigs){

        multiple_coordinates = data.frame()
        multiple_list = strsplit(x = multiple_contigs, split = "-")

        # begin for loop
        for(i in 1:length(multiple_list)){

                double_contigs = unlist(strsplit(x = multiple_list[[i]][1], split = ","))
                multiple_coordinates = rbind(multiple_coordinates, getDoubleCoordinates(chr, general, double_contigs, multiple_list[[i]][1]))

                # get the coordinates of the inverted contigs
                double_coordinates = general[general[,1] == chr & general[,4] %in% double_contigs,]
                # invert the coordinates
                double_coordinates = invertCoordinates(double_coordinates)

                multiple_coordinates = rbind(multiple_coordinates, double_coordinates[double_coordinates[,4] %in% multiple_list[[i]][2:length(multiple_list[[i]])],])

        }
        # end for loop

        # check order of coordinates
        boo = multiple_coordinates[,2] > multiple_coordinates[,3]
        multiple_coordinates[boo,] = multiple_coordinates[boo, c(1,3,2,4)]

        return(multiple_coordinates)

}
