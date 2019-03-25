def isEmpty(list):
    if not list:
        return True
    else:
        return False

def print_list(list):
	for i in range(len(list)):
		print(list[i]);

def f_comma(my_str, group, char):
    my_str = str(my_str)
    return char.join(my_str[i:i+group] for i in range(0, len(my_str), group))

def index_of_xters_in_string(string, xter):
	return [match.start() for match in re.finditer(re.escape(xter), string)]

def insert(string, array,xter):
  sequence= [0];
  new_string = "";
  i = 1;
  for i in range(len(array)):
    sequence.append(array[i]+1) #seems the regphos position annotation counts from zero
    new_string += string[sequence[i]:sequence[i+1]]+xter
  new_string += string[sequence[-1]:]+"\n"
  return new_string