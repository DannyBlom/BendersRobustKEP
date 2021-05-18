# This file is meant to convert benchmark instances of Carvalho et al. (2020) to a suitable format
function read_klimentova_graph(
        compat::String,       # file encoding compatabilities
        pairs::String,        # file encoding pairs
        output::String        # file to store the results
    )

    #= Each instances generated for the Carvalho et al. (2020) paper is determined by two files
     = 1. A file "a_b_compat.txt", where a is the number of nodes and b the instance ID
     = 2. A file "a_b_pairs.txt"
     =#
    f = open(compat)
    g = open(pairs)

    # We transform their data to a single-file format, named "Klimentova_a_b.txt"
    h = open(output, "w")

    # Read files as array of lines
    flines = readlines(f)
    glines = readlines(g)
    fcnt = 1
    gcnt = 1
    fcurrent_line = split(flines[fcnt], "\t")
    gcurrent_line = split(glines[gcnt], "\t")

    # Find number of nodes, number of arcs (notice that this includes artificial arcs from pairs to NDDs)
    nnodes = parse(Int64, fcurrent_line[1])
    narcs = parse(Int64, fcurrent_line[2])

    npairs = parse(Int64, gcurrent_line[1])
    ndonors = parse(Int64, gcurrent_line[2])

    #Check which of the nodes are donors
    donors = Set{Int64}()
    for gcnt = 2:(nnodes + 1)
        gcurrent_line = split(glines[gcnt], "\t")
        if gcurrent_line[3] == "-"
            push!(donors, parse(Int64, gcurrent_line[1]) )
        end
    end

    # Relabel nodes, such that donors have highest labels
    labels = Array{Tuple{Int64, Int64}}(undef, nnodes) #Key: label, value: corresponding node
    label = 0
    labelled_donors = 0
    for node in 0:(nnodes-1)
        if node in donors
            labels[node+1] = (node, npairs+labelled_donors)
            labelled_donors += 1
        else
            labels[node+1] = (node, label)
            label += 1
        end
    end

    # Sort back based on original nodes to make it easier to find labels
    sort!(labels, by = v -> v[1])

    # Edge information written to suitable format
    arc_list = Array{Tuple{Int64, Int64}}(undef, narcs)
    arc_ids = falses(narcs)
    for fcnt = 2:(narcs + 1)
        line = split(flines[fcnt], "\t") #Two integers indicating the head and tail of the arc
        first = parse(Int64, line[1])
        second = parse(Int64, line[2])

        if second in donors
            # artificial arcs from pairs to NDDs are omitted
            continue
        else
            label_first = labels[first+1][2]
            label_second = labels[second+1][2]
            arc_list[fcnt-1] = (label_first, label_second)
            arc_ids[fcnt-1] = true
        end
    end
    arc_list = arc_list[arc_ids]
    sort!(arc_list, by = v -> v[1])

    # write statistics to output file
    narcs = length(arc_list)
    write(h, "Nr_Pairs = $npairs\n")
    write(h, "Nr_NDD = $ndonors\n")
    write(h, "Nr_Arcs = $narcs\n")

    sort!(labels, by = v -> v[2])
    for k in 1:length(labels)
        node, node_label = labels[k]
        gcurrent_line = split(glines[node+2], "\t")

        # Write new label together with vertex failure probability, which we set to 0 for donor nodes
        if gcurrent_line[4] != "-"
            failure_prob = parse(Int64, gcurrent_line[4])/100
            write(h, "$node_label\t$(failure_prob)\n")
        else
            write(h, "$node_label\t0.0\n")
        end
    end

    for arc in arc_list
        write(h, "($(arc[1]),$(arc[2])), 0, 1\n")
    end
    close(f)
    close(g)
    close(h)
end

# main function
function  main(args)

    if length(ARGS) != 3
        print("The script expects three arguments: compatibility file, pairs file, and output file\n")
    else
        read_klimentova_graph(args[1], args[2], args[3])
    end
end

main(ARGS)
