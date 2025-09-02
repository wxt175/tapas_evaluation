###Utility Functions


#1. GPTcelltype Function
###gpt function###
gptcelltypeanno <- function(input, tissuename=NULL, model, topgenenumber = 10) {
  API.flag <- 1
  if (class(input)=='list') {
    input <- sapply(input,paste,collapse=',')
  } else {
    input <- input[input$avg_log2FC > 0,,drop=FALSE]
    input <- tapply(input$gene,list(input$cluster),function(i) paste0(i[1:topgenenumber],collapse=','))
  }
  
  if (!API.flag){
    message = paste0('Identify cell types of ',tissuename,' cells using the following markers separately for each\n row. Only provide the cell type name. Do not show numbers before the name.\n Some can be a mixture of multiple cell types. ',  "\n", paste0(names(input), ':',unlist(input),collapse = "\n"))
    
    return(message)
    
  } else {
    print("Note: OpenAI API key found: returning the cell type annotations.")
    cutnum <- ceiling(length(input)/30)
    if (cutnum > 1) {
      cid <- as.numeric(cut(1:length(input),cutnum))	
    } else {
      cid <- rep(1,length(input))
    }
    
    allres <- sapply(1:cutnum,function(i) {
      id <- which(cid==i)
      flag <- 0
      while (flag == 0) {
        k <- openai::create_chat_completion(
          model = model,
          message = list(list("role" = "user", "content" = paste0('Identify cell types of ',tissuename,' cells using the following markers separately for each\n row. Only provide the cell type name. Do not show numbers before the name.\n Some can be a mixture of multiple cell types.\n',paste(input[id],collapse = '\n'))))
        )
        res <- strsplit(k$choices[,'message.content'],'\n')[[1]]
        if (length(res)==length(id))
          flag <- 1
         Sys.sleep(2)
      }
      names(res) <- names(input)[id]
      res
    },simplify = F) 
    print('Note: It is always recommended to check the results returned by GPT in case of\n AI hallucination, before going to down-stream analysis.')
    return(gsub(',$','',unlist(allres)))
  }
  
}
###A wrapper function
CT_GPTpredict<-function(N,marker,tissueName=NULL,model){
  set.seed(2024123456)
  # Initialize a list to store results
  results_list <- list()
  
  # Run the function multiple times
  for (i in 1:N) {
    result <- gptcelltypeanno(input = marker,model = model,tissuename = tissueName
    )
    results_list[[i]] <- result
  }
  
  # Combine results column-wise using cbind
  combined_results <- do.call(cbind, results_list)
  Anno_results<-data.frame(combined_results)
  return(Anno_results)
}


#2. Function for using Chatgpt to test the relationship between two cell types, using 'Belong to' promopt
BelongTo<- function(celltype1, celltype2, model = 'gpt-4') {
  # Retrieve the OpenAI API key from environment variables
  OPENAI_API_KEY <- Sys.getenv("OPENAI_API_KEY")
  
  # Check if the API key is available
  if (OPENAI_API_KEY == "") {
    message("Note: OpenAI API key not found: returning the prompt itself.")
    API.flag <- FALSE
  } else {
    API.flag <- TRUE
  }
  
  # If API key is not found, return the generated prompt
  if (!API.flag) {
    message <- paste0('Note: OpenAI API key not found')
    
    return(message)
    
  } else {
    # If API key is available, send a request to the OpenAI API
    print("Note: OpenAI API key found: Check the entailment.")
    
    k <- openai::create_chat_completion(
      model = model,
      messages = list(
        list(
          "role" = "user", 
          "content" = paste0(
            "We are evaluating answers to the prediction of cell type name based on a list of marker genes.\n",
            "Here are two possible answers for the cell type name:\n",
            "Cell Type 1: {", celltype1, "} \n",
            "Cell Type 2: {", celltype2, "} \n",
            "Does Cell Type 1 belong to Cell Type 2? \n Reply ‘yes’ or ‘no’.")
        )
      )
    )
    
    # Extract the result from the response
    res <- strsplit(k$choices[,'message.content'],'\n')[[1]]
  }
  return(res)
}

#3. Function to classify the cell type in a hierarchical tree structure with backtracking( From top node to the bottom)
classify_cell_type_recursive <- function(new_celltype, tree_df, current_node, model = "gpt-4") {
  
  # Check if the new cell type belongs to the current node
  result <- BelongTo(new_celltype, current_node, model = model)
  
  # If the new cell type belongs to the current node, we continue checking sub-nodes
  if (result == "Yes") {
    # Find all children of the current node
    sub_nodes <- tree_df$Child[tree_df$Parent == current_node]
    
    # If there are no children, return the current node as the classification result
    if (length(sub_nodes) == 0) {
      return(current_node)  # This is a leaf node
    }
    
    # Recursively check each child node
    for (sub_node in sub_nodes) {
      classified_node <- classify_cell_type_recursive(new_celltype, tree_df, sub_node, model)
      
      # If a child node matches, return that classification
      if (!is.null(classified_node)) {
        return(classified_node)
      }
    }
    
    # If none of the children match, return the current node (backtracking)
    return(current_node)
  }
  
  # If no match found at this node, return NULL to indicate no classification at this level
  return(NULL)
}
# Wrapper function to start the classification from the top/root node
classify_cell_type_from_specific_nodes <- function(new_celltype, tree_df, starting_nodes = c("Lymphoid Lineage", "Myeloid Lineage"), model = "gpt-4") {
  
  # Iterate over the starting nodes
  for (root_node in starting_nodes) {
    classification <- classify_cell_type_recursive(new_celltype, tree_df, root_node, model)
    
    # If a classification is found, return it
    if (!is.null(classification)) {
      message(paste("The final classification for", new_celltype, "is:", classification))
      return(classification)
    }
  }
  
  # If no classification is found in any of the starting nodes, return NA
  message(paste("No match found for", new_celltype))
  return(NA)
}

# Function to classify multiple cell types from multiple lists
classify_multiple_cell_types <- function(cell_type_lists, tree_df, starting_nodes = c("Lymphoid Lineage", "Myeloid Lineage"), model = "gpt-4") {
  
  # Initialize a list to store the classification results
  classification_results <- list()
  
  # Loop over each list of cell types
  for (i in seq_along(cell_type_lists)) {
    cell_type_list <- cell_type_lists[[i]]
    
    # Loop over each cell type in the current list
    list_classification <- list()
    for (celltype in cell_type_list) {
      classification <- classify_cell_type_from_specific_nodes(celltype, tree_df, starting_nodes, model)
      list_classification[[celltype]] <- classification  # Store the classification result
    }
    
    # Add the classification results for this list to the main results list
    classification_results[[paste("List", i)]] <- list_classification
  }
  
  # Return the classification results
  return(classification_results)
}



#. 4. Assign cell type from tree by providing the tree
Assign_ct_Tree<- function(test_celltype, tree_df, model) {
  # Retrieve the OpenAI API key from environment variables
  OPENAI_API_KEY <- Sys.getenv("OPENAI_API_KEY")
  
  # Check if the API key is available
  if (OPENAI_API_KEY == "") {
    message("Note: OpenAI API key not found: returning the prompt itself.")
    API.flag <- FALSE
  } else {
    API.flag <- TRUE
  }
  
  # If API key is not found, return the generated prompt
  if (!API.flag) {
    message <- paste0('Note: OpenAI API key not found')
    
    return(message)
    
  } else {
    # Convert the tree_df to a readable string format
    tree_string <- paste(apply(tree_df, 1, function(row) {
      paste(row["Parent"], "->", row["Child"])
    }), collapse = "\n")
    # If API key is available, send a request to the OpenAI API
    print("Note: OpenAI API key found: Assign the celltype.")
    
    k <- openai::create_chat_completion(
      model = model,
      messages = list(
        list(
          "role" = "user", 
          "content" = paste0(
            "I have a predicted cell type: {", test_celltype, "}.\n",
            "Here is a tree of cell types showing parent-child relationships:\n",
            tree_string, "\n",
            "Please follow the tree structure, start from its root node, compare and assign the predicted cell type to the closest matching cell type in the tree. \n" ,
            "Return ONLY one nearest assigned cell type name, NO sentence.\n",
            "If the predicted cell type does not belong to the tree, return 'Unknown'.")
        )
      )
    )
    
    # Extract the result from the response
    res <- strsplit(k$choices[,'message.content'],'\n')[[1]]
  }
  return(res)
}
#example: assign_res<-Assign_ct_Tree("Adipocytes",pbmc_df)
# Function to classify multiple cell types from multiple lists
Assign_Multiple_Cell_Types <- function(cell_type_lists, tree_df, model) {
  # Initialize a list to store the classification results
  classification_results <- list()
  
  # Loop over each list of cell types
  for (i in seq_along(cell_type_lists)) {
    cell_type_list <- cell_type_lists[[i]]
    
    # Initialize a sub-list to store the results for this specific list
    sub_list_results <- list()
    
    # Loop over each cell type in the current list
    for (celltype in names(cell_type_list)) {
      # Assign the cell type using the Assign_ct_Tree function
      assigned_celltype <- Assign_ct_Tree(celltype, tree_df, model)
      
      # Store the result in the sub-list
      sub_list_results[[celltype]] <- assigned_celltype
      # Add a delay to avoid hitting the rate limit
      Sys.sleep(3)  # Wait for 2 seconds between requests (adjust as necessary)
    }
    
    # Add the results of the current list to the main results list
    classification_results[[paste("List", i)]] <- sub_list_results
  }
  
  # Return the overall classification results
  return(classification_results)
}


#5. Process the assigned results list
Process_AssignRes <- function(assigned_result_list,frequncy_table) {
  # Initialize a list to store the processed data frames
  processed_results <- list()
  for (i in seq_along(assigned_result_list)) {
    assigned_celltypes <- assigned_result_list[[i]]
    assigned_celltypes_vector <- unlist(assigned_celltypes)
    freq_table<-frequncy_table[[i]]
    freq_vector <- unlist(freq_table)
    freq_vector<-as.numeric(freq_vector)
    frequency_table <- data.frame(cbind(assigned_celltypes_vector,freq_vector))
    
    colnames(frequency_table) <- c("Assigned_Cell_Type", "Frequency")
    frequency_table$Frequency<-as.numeric(frequency_table$Frequency)
    node_count <- aggregate(Frequency ~ Assigned_Cell_Type, data = frequency_table, sum)
    # Store the processed table for this list
    processed_results[[paste("List", i)]] <- node_count
  }
  
  return(processed_results)
}

### 6. Get the distance matrix for the assigned node
Get_matrix<-function(assigned_result_list,distance_mtx,idx){
  Node_distance <- list()
  for (i in idx) {
  assigned_celltypes <- assigned_result_list[[i]]
  celltype<- unique(unlist(assigned_celltypes))
  if(length(celltype)==1){
    sub_dis_mtx=0
  }else{
    sub_dis_mtx <- distance_mtx[celltype,celltype]
   }
  Node_distance[[i]] <- sub_dis_mtx
 }
  return(Node_distance)
}


###7. Calculate the average weight matrix
average_distance <- function(distance_matrix, n_sample, N = 20,D_Unknown) {
  # Extract specific nodes based on the provided sample indices
  specific_nodes <- rownames(distance_matrix)[n_sample]
  repetition_df <- as.data.frame(table(specific_nodes))
  colnames(repetition_df) <- c('cell_type', 'Freq')
  unique_node <- repetition_df$cell_type
  
  # Initialize variable to store the average distance
  weighted_average_distance <- NA
  
  # Handle cases with only one unique node
  if (length(unique_node) == 1) {
    if (unique_node == 'Unknown') {
      weighted_average_distance <- D_Unknown  # Special case for 'Unknown'
    } else {
      weighted_average_distance <- 0  # If only one node and not 'Unknown'
    }
  } else {
    # Subset the distance matrix for the specific nodes, excluding "Unknown" if it exists
    filtered_nodes <- unique_node[unique_node != "Unknown"]
    
    if (length(filtered_nodes) > 0) {
      specific_distance_matrix <- distance_matrix[filtered_nodes, filtered_nodes, drop = FALSE]
      
      # Create a named vector for repetition counts
      specific_repetition_counts <- setNames(repetition_df$Freq, repetition_df$cell_type)
      
      # Initialize total weighted distance and total weight
      total_weighted_distance <- 0
      total_weight <- 0.5 * (N - 1) * N  # Based on the formula
      
      # Calculate the weighted distance for known nodes
      for (i in 1:nrow(specific_distance_matrix)) {
        for (j in i:ncol(specific_distance_matrix)) {
          # Get the pairwise distance between nodes i and j
          distance <- specific_distance_matrix[i, j]
          
          # Get the repetition counts for nodes i and j using node names
          rep_i <- specific_repetition_counts[filtered_nodes[i]]
          rep_j <- specific_repetition_counts[filtered_nodes[j]]
          
          # Calculate the weighted distance
          weighted_distance <- distance * rep_i * rep_j
          
          # Add the weighted distance to the total
          total_weighted_distance <- total_weighted_distance + weighted_distance
        }
      }
      
      # Calculate the weighted distance for "Unknown" nodes if they exist
      n_unknown <- specific_repetition_counts["Unknown"]
      D_unknown <- D_Unknown
      if (!is.na(n_unknown) && n_unknown > 0) {
        # Add the distance for "Unknown" nodes (self-distances)
        total_distance <- total_weighted_distance + (D_unknown * 0.5 * n_unknown * (n_unknown - 1))
        
        # Add distances between "Unknown" nodes and all other nodes
        for (i in 1:length(filtered_nodes)) {
          rep_i <- specific_repetition_counts[filtered_nodes[i]]
          total_distance <- total_distance + (D_unknown * rep_i * n_unknown)
        }
      } else {
        total_distance <- total_weighted_distance
      }
      
      # Calculate the weighted average distance
      weighted_average_distance <- as.numeric(total_distance / total_weight)
    }
  }
  
  return(weighted_average_distance)  # Return the average distance
}





Get_avg_dis <- function(assigned_result_list, frequency_table_list, distance_mtx, N = 20,D_Unknown) {
  Avg_distance <- list()  # Initialize list to store average distances
  
  # Loop through each assigned result
  for (idx in seq_along(assigned_result_list)) {
    assigned_celltypes <- assigned_result_list[[idx]]
    assigned_celltypes_vector <- unlist(assigned_celltypes)
    
    freq_table <- frequency_table_list[[idx]]
    freq_vector <- as.numeric(unlist(freq_table))
    
    # Create a frequency table
    frequency_table <- data.frame(cbind(assigned_celltypes_vector, freq_vector))
    colnames(frequency_table) <- c("Assigned_Cell_Type", "Frequency")
    frequency_table$Frequency <- as.numeric(frequency_table$Frequency)
    
    # Aggregate frequencies to get node counts
    node_count <- aggregate(Frequency ~ Assigned_Cell_Type, data = frequency_table, sum)
    
    # Extract the specific nodes
    specific_nodes <- node_count$Assigned_Cell_Type
    
    # Check if all specific nodes are present in the distance matrix
    if (!all(specific_nodes %in% rownames(distance_mtx)) || !all(specific_nodes %in% colnames(distance_mtx))) {
      # If any specific node is not found, return NA
      weighted_average_distance <- NA
    } else if (length(specific_nodes) == 1) {
      # Handle cases with only one specific node
      if (specific_nodes == 'Unknown') {
        weighted_average_distance <- D_Unknown  # Special case for 'Unknown'
      } else {
        weighted_average_distance <- 0  # If only one node and not 'Unknown'
      }
    } else {
      # Subset the distance matrix for the specific nodes, excluding "Unknown" if it exists
      filtered_nodes <- specific_nodes[specific_nodes != "Unknown"]
      specific_distance_matrix <- distance_mtx[filtered_nodes, filtered_nodes, drop = FALSE]
      
      # Create a named vector for repetition counts
      specific_repetition_counts <- setNames(node_count$Frequency, node_count$Assigned_Cell_Type)
      
      # Initialize total weighted distance and total weight
      total_weighted_distance <- 0
      total_weight <- 0.5 * (N - 1) * N  # Based on the formula
      
      # Calculate the weighted distance for known nodes
      for (i in 1:nrow(specific_distance_matrix)) {
        for (j in i:ncol(specific_distance_matrix)) {
          # Get the pairwise distance between nodes i and j
          distance <- specific_distance_matrix[i, j]
          
          # Get the repetition counts for nodes i and j using node names
          rep_i <- specific_repetition_counts[filtered_nodes[i]]
          rep_j <- specific_repetition_counts[filtered_nodes[j]]
          
          # Calculate the weighted distance
          weighted_distance <- distance * rep_i * rep_j
          
          # Add the weighted distance to the total
          total_weighted_distance <- total_weighted_distance + weighted_distance
        }
      }
      
      # Calculate the weighted distance for "Unknown" nodes if they exist
      n_unknown <- specific_repetition_counts["Unknown"]
      D_unknown <- D_Unknown
      if (!is.na(n_unknown) && n_unknown > 0) {
        # Add the distance for "Unknown" nodes (self-distances)
        total_distance <- total_weighted_distance + (D_unknown * 0.5 * n_unknown * (n_unknown - 1))
        
        # Add distances between "Unknown" nodes and all other nodes
        for (i in 1:length(filtered_nodes)) {
          rep_i <- specific_repetition_counts[filtered_nodes[i]]
          total_distance <- total_distance + (D_unknown * rep_i * n_unknown)
        }
      } else {
        total_distance <- total_weighted_distance
      }
      
      # Calculate the weighted average distance
      weighted_average_distance <- total_distance / total_weight
    }
    
    # Store the result for this iteration
    Avg_distance[[idx]] <- as.numeric(weighted_average_distance)
  }
  
  return(Avg_distance)  # Return the list of average distances
}



###calculate the p-value
get_pval<-function(score,shape,rate){
  # Lower tail (P(X <= x)) for values less than or equal to x
  p_value_lower <- pgamma(score, shape = shape, rate = rate)
  # Upper tail (P(X >= x)) for values greater than or equal to x
  p_value_upper <- 1 - p_value_lower
  return(p_value_upper)
}

get_p_vector<-function(vector,shape,rate){
  p_res<-numeric()
  for (i in 1:length(vector)){
    p_res[i]<-get_pval(vector[i],shape = shape, rate = rate)
  }
  return(p_res)
}



