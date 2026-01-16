const global where_my_CPTS_are = @__DIR__

function getMyCPTsCollectionPath()
    return where_my_CPTS_are
end

function getMyCPTsList()
    return readdir(getMyCPTsCollectionPath())
end

function getMyCPTPath(cpt_name)
    # Get the path only if the CPT exists
    if cpt_name in getMyCPTsList()
         raw_path = joinpath(getMyCPTsCollectionPath(), cpt_name)
        return "\"$raw_path\""  # Handles possible spaces in path
    else
        print("Your cpt does not exist in path!")
        return nothing
    end
end

function lookCPTcontent(cpt_name)
    if cpt_name in getMyCPTsList()
        raw_path = joinpath(getMyCPTsCollectionPath(), cpt_name)
        return readlines(raw_path)
    else
        println("Your cpt does not exist in path!")
        return nothing
    end
end