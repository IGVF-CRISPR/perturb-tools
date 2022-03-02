def _annotate_homopolymer(df):
    
    """Annotate any homopolymer sequences in sgRNAs"""
    
    df['polyT'] = df.protospacer.str.contains("TTTT")
    df['polyG'] = df.protospacer.str.contains("GGGG")
    df['polyC'] = df.protospacer.str.contains("CCCC")
    df['polyA'] = df.protospacer.str.contains("AAAA")
    
    return df

def _get_GC_content(df):
    
    """Get the GC content of sgRNAs"""
    
    g_count = df.protospacer.str.count("G")
    c_count = df.protospacer.str.count("C")
    
    df['GC_content'] = (g_count + c_count)*100/20
    
    return df

def _annotate_guide_biochemistry(df):
    
    """Adds information to the guide df about GC content and polyT content"""
    
    df = _annotate_homopolymer(df)
    df = _get_GC_content(df)
    
    return df