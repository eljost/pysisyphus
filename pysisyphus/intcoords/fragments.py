def merge_fragments(fragments):
    """Merge a list of sets."""
    # Hold the final fragments that can't be merged further, as they
    # contain distinct atoms.
    merged = list()
    while len(fragments) > 0:
        popped = fragments.pop(0)
        # Look for an intersection between the popped unmerged fragment
        # and the remaining unmerged fragments.
        for frag in fragments:
            if popped & frag:
                fragments.remove(frag)
                # If a intersecting unmerged fragment is found merge
                # both fragments and append them at the end.
                fragments.append(popped | frag)
                break
        else:
            # Add the unmerged fragment into merged if it doesn't
            # intersect with any other unmerged fragment.
            merged.append(popped)
    return merged
