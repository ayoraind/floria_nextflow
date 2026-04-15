// include { check_mandatory_parameter; check_optional_parameters; check_parameter_value } from './params_utilities.nf'

def default_params(){
    /***************** Setup inputs and channels ************************/
    def params = [:]
    // Defaults for configurable variables
    params.help = false
    params.version = false
   //  params.publish_dir_mode = 'copy'
    params.reads = false
    params.assemblies = false
    params.output_dir = false
    return params
}

def check_params(Map params) { 
    def final_params = params
    
    // set up reads files
    final_params.reads = check_mandatory_parameter(params, 'reads')
     
    // set up output directory
    final_params.output_dir = check_mandatory_parameter(params, 'output_dir') - ~/\/$/
     
    // set up assemblies filepath
    final_params.assemblies = check_mandatory_parameter(params, 'assemblies')
    	
    return final_params
}

def check_mandatory_parameter(Map params, String parameter_name){
    if ( !params[parameter_name]){
        println "You must specify a " + parameter_name
        System.exit(1)
    } else {
        return params[parameter_name]
    }
}

def check_optional_parameters(Map params, List parameter_names){
    if (parameter_names.collect{name -> params[name]}.every{param_value -> param_value == false}){
        println "You must specifiy at least one of these options: " + parameter_names.join(", ")
        System.exit(1)
    }
}

def check_parameter_value(String parameter_name, String value, List value_options){
    if (value_options.any{ it == value }){
        return value
    } else {
        println "The value for " + parameter_name + " must be one of " + value_options.join(", ")
        System.exit(1)
    }
}

def rename_params_keys(Map params_to_rename, Map old_and_new_names) {
    old_and_new_names.each{ old_name, new_name ->
        if (params_to_rename.containsKey(old_name))  {
            params_to_rename.put( new_name, params_to_rename.remove(old_name ) )
        }
    }
    return params_to_rename
}
