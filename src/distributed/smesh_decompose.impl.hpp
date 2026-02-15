// Decomposition algorithm for distributed mesh:

// 1) Given the communicator size P determin ownership ranges of elements so that each rank has (Ne / P), the reminder of Ne / P goes to the last rank
// 2) For the nodes we need to have a temporary decomposition Ln = Nn / P ...
// 3) Determine how many elements have nodes in each partition count elements[d][i] / Ln
// 4) Create node to element index (n2e) for each parition, we need, options:
//      a. Unique node local 2 global node indexing by counting and bucketing (owned + ghosts) ?
//      b. Element to node (partials) to be sent
// 5) Send the partial n2e to the node's owner rank (3 arrays: node_idx, n2e_ptr, n2e_idx)
// 6) Compute the full n2e in the temp owner rank
// 7) Elements are replicated for each connection, we can now:
//      a. Define node actual ownership and determine ghosts, owner of the node is the smallest rank  
//         with  an incident element
//      b. Identify aura elements by using the full n2e graph
//      c. Compute partial n2n graph
// 8) Send full n2e to the ranks owning the elements (node_idx, n2e_ptr, n2e_idx), this  gives:
//      a. The full n2e graph for each rank with incident elements
// 9) Send n2n partials to temp owner rank and compute the full n2n graph
// 10) Having the n2e we can send the n2n parts to the ranks with incident elements


// For option 4.a
// Counting/Bucketing algorithm for node_mapping (local2global)
// 1) Determine ranges of nodes [begin, end) for each partition
// 2) For each non-empy partition create buckets
// 3) For each bucket determine range and use counting sort to create an uniquely sorted list of nodes
// 4) Create partial n2e


// When I have the n2e we can: 
// - for each node get the element and replace the global index with the local index (no sorting required)