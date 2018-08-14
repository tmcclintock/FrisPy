#This function parses a line in the configuration file
#It pulls out a key and a value and adds them to the dictionary
#that is config.
def parse_variable(config, var):
    #Make sure the line is formatted properly
    if '=' not in var:
        message = "Improper variable: %s. Use syntax: key = value."%var
        raise Exception(message)
    key, value = var.split('=',1)
    key = key.strip() #Cuts off whitespace
    #Cut off a trailing comment
    if '#' in value: value = value.split('#')[0]
    value = value.strip()
    #Add the key and value pair
    config[key] = value

#Construct the config dictionary with key/value pairs
def read_config_file(file_name):
    config = dict()
    file_in = open(file_name)
    for line in file_in:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            pass
        else:
            parse_variable(config,line)

    file_in.close()
    return config
