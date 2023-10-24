#!/usr/bin/env python3
# Wrap all the F77 functions of the commprof.cpp with extern "C"

from itertools import takewhile
#def wrap_functions_with_extern_c(content):
#    pattern = r"(void F77[A-Z0-9_]+\(.+?\)\n\{.+?\n\})"
#    def wrap_with_extern_c(match):
#        return 'extern "C" {\n' + match.group(1) + "\n}"
#    wrapped_content = re.sub(pattern, wrap_with_extern_c, content, flags=re.DOTALL)
#    return wrapped_content
def wrap_functions_with_extern_c(content):
    lines = content.split('\n')
    new_content = []
    i = 0

    while i < len(lines):
        # If current line is "void" and the next line starts with "F77", then it's a potential target
        if lines[i].strip() == "void" and i + 1 < len(lines) and lines[i + 1].strip().startswith("F77"):
            # Check if it's already wrapped with extern "C"
            if i > 0 and lines[i - 1].strip() == 'extern "C" {':
                new_content.append(lines[i])
            else:
                # Add the extern "C" line before the function
                indentation = ''.join(list(takewhile(str.isspace, lines[i])))
                new_content.append(indentation + 'extern "C" {')
                new_content.append(lines[i])

                # Increment i and add all the function lines until its closing bracket
                i += 1
                while i < len(lines) and not lines[i].strip().startswith("}"):
                    new_content.append(lines[i])
                    i += 1

                # Append the closing bracket of the function and '}' to close the extern "C"
                if i < len(lines):
                    new_content.append(lines[i])
                new_content.append(indentation + "}")
        else:
            new_content.append(lines[i])
        i += 1

    return '\n'.join(new_content)


if __name__ == '__main__':
    with open('old_commprof.cpp', 'r') as file:
        content = file.read()

    # Get the wrapped content
    new_content = wrap_functions_with_extern_c(content)

    # Write to a new file
    with open('commprof.cpp', 'w') as file:
        file.write(new_content)

