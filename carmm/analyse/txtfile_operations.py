LINE_NOT_FOUND = object()

def reverse_search_for(keys, lines, line_start=0):
    for ll, line in enumerate(lines[line_start:][::-1]):
        if any(key in line for key in keys):
            return len(lines) - ll - 1

    return LINE_NOT_FOUND


