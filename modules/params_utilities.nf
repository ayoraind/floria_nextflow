include { version_message; help_message } from './messages.nf'

def help_or_version(Map params, String version){
    // Show help message
    if (params.help){
        version_message(version)
        help_message()
        System.exit(0)
    }

    // Show version number
    if (params.version){
        version_message(version)
        System.exit(0)
    }
}

