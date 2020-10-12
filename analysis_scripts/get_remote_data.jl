function remote_download(server, directory, file)
    commandline = server * directory * file
    run(`rsync -av $commandline .`)
    return nothing
end