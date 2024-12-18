import argparse
from pathlib import Path
from urllib.parse import urljoin, urlsplit
import time
from tqdm import tqdm

import requests
from requests.exceptions import JSONDecodeError as RequestsJSONDecodeError

from requests_toolbelt.multipart.encoder import (
    MultipartEncoder,
    MultipartEncoderMonitor,
)



from xchemalign import utils


LOGIN_URL = "/accounts/login/"
UPLOAD_URL = "/api/upload_target_experiments/"
LANDING_PAGE_URL = '/viewer/react/landing/'


# this needs to be kept more or less up to date
USER_AGENT = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"


logger = utils.Logger()


def check_deployment(session, auth_token):

    session.get(url)  # sets cookie

    csrftoken = session.cookies.get('csrftoken', '')
    session.headers.update({
        "X-CSRFToken": csrftoken,
        "User-Agent": USER_AGENT,
    })

    if auth_token:
        session.cookies.update({
            "sessionid": "",
        })
     
    
    return session
    
    
    


def upload(input_file, url, proposal, auth_token=None):
    splits = urlsplit(url)
    base_url = f'{splits.scheme}://{splits.netloc}'
    upload_url = urljoin(base_url, UPLOAD_URL)

    filename = Path(input_file).name
    
    with requests.Session() as session:
        
        session.headers.update({
            "User-Agent": USER_AGENT,
            "Referer": urljoin(url, LANDING_PAGE_URL),
            "Referrer-policy": "same-origin",
        })          

        # session = check_deployment(session, auth_token)
        session.get(url)  # sets csrftoken
    
        csrftoken = session.cookies.get('csrftoken', None)
        if csrftoken:
            session.headers.update({
                "X-CSRFToken": csrftoken,
                "User-Agent": USER_AGENT,
            })
    
        if auth_token:
            session.cookies.update({
                "sessionid": "tui2tg6ic6kkmcwzzo734lz10oad21zq",
            })        
       
        encoder = MultipartEncoder(
            fields={
                "target_access_string": proposal,
                "file": (
                    input_file,
                    open(input_file, "rb"),
                    "application/octet-stream",
                ),
            }
        )
       
        with open(input_file, "rb") as f:
            file_size = int(f.seek(0, 2)) 
            f.seek(0)

            def callback(monitor):
                pbar.update(monitor.bytes_read - pbar.n)

            monitor = MultipartEncoderMonitor(encoder, callback)

            session.headers.update({"Content-Type": monitor.content_type})

            with tqdm(
                total=file_size, unit="B", unit_scale=True, desc='Uploading'
            ) as pbar:
                response = session.post(
                    upload_url,
                    data=monitor,
                )            

        if response.ok:
            # not quite there yet
            try:
                jresp = response.json()
                # and once more
                try:
                    task_status_url = jresp["task_status_url"]
                except KeyError:
                    logger.error('key error, what does it mean?')
                    return
            except requests.exceptions.JSONDecodeError:
                logger.error('json error what does it mean?')
                return

            # should get back "task_status_url" in response. keep polling for updates
            task_status = {}
            count = 0
            status_url = urljoin(base_url, task_status_url)
            finished = False
            errors = []
            
            logger.info(f'{filename} uploaded, processing:')
            while True:
                # ping the task url endpoint for a while and print the messages
                status = session.get(status_url)
                count = count + 1
                try:
                    statjson = status.json()
                except RequestsJSONDecodeError:
                    continue
                
                if task_status != statjson:
                    task_status = statjson
                    msg = task_status.get("messages", None)
                    if isinstance(msg, list):
                        for m in msg:
                            logger.info(m)
                            if m.startswith('ERROR'):
                                errors.append(m)
                    else:
                        logger.info(f'{task_status.get("status", "")}: {msg}')
                    if task_status in ("SUCCESS", "FAILED", "CANCELED", "FATAL"):
                        finished = True
                        break
                    
                if statjson.get("ready", False) == True:
                    # if task complete, quit
                    finished = True
                    break
                if statjson.get("status", "") in ("SUCCESS", "FAILED", "CANCELED", "FATAL"):
                    # if task complete, quit
                    finished = True
                    break                    
                if count > 100:
                    # if no changes for {randomly selected number of pings}, quit
                    break

                time.sleep(1)

            if finished:
                if errors:
                    logger.info(f'The following errors were encountered when processing {filename}')
                    for e in errors:
                        logger.info(e)
                
        else:
            # upload response else, failure
            print("Upload failed: error {}, {}".format(response.status_code, response.reason))  


        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Upload files to fragalysis")
    parser.add_argument("-i", "--input", required=True, help="Input")
    parser.add_argument("-u", "--url", required=True, help="Url to upload")
    parser.add_argument("-p", "--proposal", required=True, help="Project/proposal")
    parser.add_argument("-t", "--token", default=None, required=False, help="Authentication token")

    args = parser.parse_args()

    upload(args.input, args.url, args.proposal, auth_token=args.token)

    # normal upload
    # python src/xchemalign/uploader.py -i ../test_data/A71EV2A_xca_staging_20241104_fake_aliases.tar.gz -u http://localhost:8080/ -p lb18145-1    
    
