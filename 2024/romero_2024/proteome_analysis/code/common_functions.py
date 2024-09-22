def split_list_equally(input_list, n):
    """ Split input_list into n equally-sized sublists. """
    # Calculate the size of each sublist
    size = len(input_list) // n
    # Calculate any remainder to distribute
    remainder = len(input_list) % n
    
    # Initialize starting index for slicing
    start = 0
    result = []
    
    for i in range(n):
        # Calculate slice end index
        end = start + size + (1 if i < remainder else 0)
        # Append sublist to result
        result.append(input_list[start:end])
        # Update start index for next sublist
        start = end
    
    return result

#uid2s_over_t ends up as a dictionary where the order is defined by the S_over_T mean value smallest to largest
def annotate_with_quantiles(yp, q=8, proteins_with_domains=True):
    """
    Takes in a yeast proteome object and annotates each protein with a quantile value.

    If proteins_with_domains is True, only proteins with at least one domain are considered.

    Returns a list of quantile values, one for each protein in the proteome, as well as
    annotating each protein with a 'quantile' attribute which is equal to one of the
    quantile values.

    Expectation is that EVERY protein passed has a S_over_T_mean attribute at the 
    protein level.

    Parameters
    ----------
    yp : Proteome
        A SHEPHARD proteome object.

    q : int
        The number of quantiles to split the proteins into. Default is 8.

    proteins_with_domains : bool
        If True, only proteins with at least one domain are considered. 
        Default is True.

    Returns
    -------
    list
        A list of quantile values; every protein passed will be annotated with a quantile 
        value as an attribute that matches one of these quantile values.

    """

    # initialize set of protein UIDs
    hits = set([])

    # if proteins_with_domains is true, we preselect only for proteins with domains
    if proteins_with_domains:

        # iterate over all domains in the proteome
        for d in yp.domains:

            # extract the protein UID and add it to the set
            uid = d.protein.unique_ID
            if uid not in hits:
                hits.add(uid)
    # if proteins_with_domains is false, we just consider every protein regardless
    # of whether it has a domain or not
    else:
        for p in yp:
            hits.add(p.unique_ID)

    # create a dictionary of protein UIDs and their S_over_T_mean values
    uid2s_over_t = {}
    for h in hits:
        uid2s_over_t[h] = yp.protein(h).attribute('S_over_T_mean')

    # sort the dictionary by S_over_T_mean with smallest to largest (least soluble to most soluble)
    uid2s_over_t = dict(sorted(uid2s_over_t.items(), key=lambda item: item[1], reverse=False))
    
    # get lists of protein UIDs ordered by S_over_T_mean
    ordered_uids = list(uid2s_over_t.keys())
    
    # split into q different equally-sized groups
    B = split_list_equally(ordered_uids,8)

    # sanity check - print the number of proteins in each group to the scnree
    print('Sanity check number of proteins in each group')
    print([len(x) for x in B])

    # define the quantile value for the first group as 0.5* first quantile
    quantile = 0 + int(100*(1/q)/2)
    all_quantiles = []

    # iterate over each group of proteins
    for sublist in B:

        # iterate over each protein in the group
        all_quantiles.append(quantile)
        for uid in sublist:
            p = yp.protein(uid)
            p.add_attribute('quantile', quantile, safe=False)

        # increment the quantile value
        quantile = quantile + int(100*(1/q))

    return all_quantiles
