function downloadHCPData

datafile = './data/hcp_rsfMRI.mat';

if ~exist(datafile,'file')
  article_url = '<a href = "https://figshare.com/articles/hcp_rsfMRI_mat/12084636">https://figshare.com/articles/hcp_rsfMRI_mat/12084636</a>';
  download_url = 'https://ndownloader.figshare.com/files/22217589';
  
  fprintf('HCP data not found in the expected location: %s.\n', datafile);
  
  fprintf('You can download it here: %s\n', article_url);
  prompt = 'Do you want me to download it now? Y/N [Y]: ';
  reply = input(prompt,'s');
  
  if isempty(reply)
    reply = 'Y';
  end
  
  if upper(reply) ~= 'Y'
    fprintf('OK, exiting program now.\n');
    return
  end
  
  if ~exist('data','dir')
    mkdir('data');
  end
  
  fprintf('OK, getting it now. This may take a while...\n');
  
  % Download the hcp_rsfMRI.mat file from figshare
  out_fname = websave(datafile,download_url);
  fprintf(['Done. Human Connectome Project rsfMRI data ',...
            'downloaded from figshare to:\n%s\n'],out_fname);
else
  fprintf('Data file %s exists.\n', datafile);
end